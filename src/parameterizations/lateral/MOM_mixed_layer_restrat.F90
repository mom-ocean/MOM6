!> This module implements a parameterization of unresolved viscous
!!  mixed layer restratification of the mixed layer as described in
!!  Fox-Kemper, Ferrari and Hallberg (JPO, 2008), and
!!  whose impacts are described in Fox-Kemper et al.
!!  (Ocean Modelling).

module MOM_mixed_layer_restrat

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,     only : hchksum
use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_diag_mediator, only : diag_update_remap_grids
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_forcing_type,  only : forcing
use MOM_grid,          only : ocean_grid_type
use MOM_hor_index,     only : hor_index_type
use MOM_io,            only : vardesc, var_desc
use MOM_restart,       only : register_restart_field, MOM_restart_CS
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_EOS,           only : calculate_density

implicit none ; private

#include <MOM_memory.h>

public mixedlayer_restrat
public mixedlayer_restrat_init
public mixedlayer_restrat_register_restarts

!> Control structure for module
type, public :: mixedlayer_restrat_CS ; private
  real    :: ml_restrat_coef       !<  A non-dimensional factor by which the
                                   !! instability is enhanced over what would be
                                   !! predicted based on the resolved  gradients.  This
                                   !! increases with grid spacing^2, up to something
                                   !! of order 500.
  logical :: MLE_use_PBL_MLD       !< If true, use the MLD provided by the PBL parameterization.
                                   !! if false, MLE will calculate a MLD based on a density difference
                                   !! based on the parameter MLE_DENSITY_DIFF.
  real    :: MLE_MLD_decay_time    !< Time-scale to use un running-mean when MLD is retreating (s).
  real    :: MLE_density_diff      !< Density difference used in detecting mixed-layer
                                   !! depth (kg/m3).
  real    :: MLE_tail_dh           !< Fraction by which to extend the mixed-layer re-stratification
                                   !! depth used for a smoother stream function at the base of
                                   !! the mixed-layer.
  logical :: debug = .false.       !< If true, calculate checksums of fields for debugging.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  real, dimension(:,:), pointer :: &
         MLD          => NULL(), & !< Mixed layer depth used in the MLE re-stratification parameterization (H units)
         MLD_filtered => NULL()    !< Time-filtered MLD (H units)

  integer :: id_urestrat_time
  integer :: id_vrestrat_time
  integer :: id_uhml = -1
  integer :: id_vhml = -1
  integer :: id_MLD  = -1
  integer :: id_Rml  = -1
  integer :: id_uDml = -1
  integer :: id_vDml = -1
  integer :: id_uml  = -1
  integer :: id_vml  = -1

end type mixedlayer_restrat_CS

character(len=40)  :: mod = "MOM_mixed_layer_restrat"  ! This module's name.

contains

!>  This subroutine does interface depth diffusion.  The fluxes are
!!  limited to give positive definiteness, and the diffusivities are
!!  limited to guarantee stability.
subroutine mixedlayer_restrat(h, uhtr, vhtr, tv, fluxes, dt, MLD, G, GV, CS)
  type(ocean_grid_type),                     intent(inout) :: G      !< ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h      !< layer thickness (H units = m or kg/m2)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uhtr   !< accumulated zonal mass flux (m3 or kg)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vhtr   !< accumulated merid mass flux (m3 or kg)
  type(thermo_var_ptrs),                     intent(in)    :: tv     !< thermodynamic variables structure
  type(forcing),                             intent(in)    :: fluxes !< pointers to forcing fields
  real,                                      intent(in)    :: dt     !< time increment (sec)
  real, dimension(:,:),                      pointer       :: MLD    !< Mixed layer depth provided by PBL (H units)
  type(mixedlayer_restrat_CS),               pointer       :: CS     !< module control structure

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "Module must be initialized before it is used.")

  if (GV%nkml>0) then
    call mixedlayer_restrat_BML(h, uhtr, vhtr, tv, fluxes, dt, G, GV, CS)
  else
    call mixedlayer_restrat_general(h, uhtr, vhtr, tv, fluxes, dt, MLD, G, GV, CS)
  endif

end subroutine mixedlayer_restrat


!> This subroutine does interface depth diffusion.  The fluxes are
!! limited to give positive definiteness, and the diffusivities are
!! limited to guarantee stability.
subroutine mixedlayer_restrat_general(h, uhtr, vhtr, tv, fluxes, dt, MLD, G, GV, CS)
  type(ocean_grid_type),                     intent(inout) :: G       !< ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV      !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h       !< layer thickness (H units = m or kg/m2)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uhtr    !< accumulated zonal mass flux (m3 or kg)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vhtr    !< accumulated merid mass flux (m3 or kg)
  type(thermo_var_ptrs),                     intent(in)    :: tv      !< thermodynamic variables structure
  type(forcing),                             intent(in)    :: fluxes  !< pointers to forcing fields
  real,                                      intent(in)    :: dt      !< time increment (sec)
  real, dimension(:,:),                      pointer       :: MLD     !< Mixed layer depth provided by PBL (H units)
  type(mixedlayer_restrat_CS),               pointer       :: CS      !< module control structure

  real :: uhml(SZIB_(G),SZJ_(G),SZK_(G)) ! zonal mixed layer transport (m3/s or kg/s)
  real :: vhml(SZI_(G),SZJB_(G),SZK_(G)) ! merid mixed layer transport (m3/s or kg/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    h_avail               ! The volume available for diffusion out of each face of each
                          ! sublayer of the mixed layer, divided by dt, in units
                          ! of H * m2 s-1 (i.e., m3 s-1 or kg s-1).
  real, dimension(SZI_(G),SZJ_(G)) :: &
    htot, &               ! The sum of the thicknesses of layers in the mixed layer (H units)
    Rml_av                ! g_Rho0 times the average mixed layer density (m s-2)
  real :: g_Rho0          ! G_Earth/Rho0 (m4 s-2 kg-1)
  real :: Rho0(SZI_(G))   ! Potential density relative to the surface (kg m-3)
  real :: p0(SZI_(G))     ! A pressure of 0 (Pa)

  real :: h_vel           ! htot interpolated onto velocity points in metre (not H).
  real :: absf            ! absolute value of f, interpolated to velocity points (s-1)
  real :: u_star          ! surface friction velocity, interpolated to velocity points (m s-1)
  real :: mom_mixrate     ! rate at which momentum is homogenized within mixed layer (s-1)
  real :: timescale       ! mixing growth timescale (sec)
  real :: h_neglect       ! tiny thickness usually lost in roundoff so can be neglected (H units)
  real :: dz_neglect      ! A tiny thickness (in m) that is usually lost in roundoff so can be neglected
  real :: I4dt            ! 1/(4 dt) (sec-1)
  real :: Ihtot           ! total mixed layer thickness
  real :: a(SZK_(G))      ! A nondimensional value relating the overall flux
                          ! magnitudes (uDml & vDml) to the realized flux in a
                          ! layer.  The vertical sum of a() through the pieces of
                          ! the mixed layer must be 0.
  real :: uDml(SZIB_(G))  ! The zonal and meridional volume fluxes in the upper
  real :: vDml(SZI_(G))   ! half of the mixed layer in H m2 s-1 (m3 s-1 or kg s-1).
  real :: utimescale_diag(SZIB_(G),SZJ_(G)) ! restratification timescales
  real :: vtimescale_diag(SZI_(G),SZJB_(G)) ! in the zonal and meridional
                                            ! directions, in s, stored in 2-D
                                            ! arrays for diagnostic purposes.
  real :: uDml_diag(SZIB_(G),SZJ_(G)), vDml_diag(SZI_(G),SZJB_(G))
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  real, dimension(SZI_(G)) :: rhoSurf, deltaRhoAtKm1, deltaRhoAtK, dK, dKm1, pRef_MLD ! Used for MLD
  real, dimension(SZI_(G)) :: rhoAtK, rho1, d1, pRef_N2 ! Used for N2
  real :: aFac, bFac, ddRho
  real :: hAtVel, zIHaboveVel, zIHbelowVel

  real :: PSI, PSI1, z ! For statement function

  ! Stream function as a function of non-dimensional position within mixed-layer (F77 statement function)
  !PSI(z) = max(0., (1. - (2.*z+1.)**2 ) )
  PSI1(z) = max(0., (1. - (2.*z+1.)**2 ) * (1. + (5./21.)*(2.*z+1.)**2) )
  real :: BOTTOP, XP, PSIC, DD
  BOTTOP(z) = 0.5*(1.-SIGN(1.,z+0.5)) ! =0 for z>-0.5, =1 for z<-0.5
  XP(z) = max(0., min(1., (-z-0.5)*2./(1.+2.*CS%MLE_tail_dh) ) )
  DD(z) = (1.-3.*(XP(z)**2)+2.*(XP(z)**3))**(1.+2.*CS%MLE_tail_dh)
  PSIC(z) = max( PSI1(z), DD(z)*BOTTOP(z) )
  PSI(z) = PSIC(z)

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not.associated(tv%eqn_of_state)) call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "An equation of state must be used with this module.")

  if (CS%MLE_density_diff > 0.) then ! We need to calculate a mixed layer depth, MLD.
    !! TODO: use derivatives and mid-MLD pressure. Currently this is sigma-0. -AJA
    pRef_MLD(:) = 0.
    do j = js-1, je+1
      dK(:) = 0.5 * h(:,j,1) * GV%H_to_m ! Depth of center of surface layer
      call calculate_density(tv%T(:,j,1), tv%S(:,j,1), pRef_MLD, rhoSurf, is-1, ie-is+3, tv%eqn_of_state)
      deltaRhoAtK(:) = 0.
      CS%MLD(:,j) = 0.
      do k = 2, nz
        dKm1(:) = dK(:) ! Depth of center of layer K-1
        dK(:) = dK(:) + 0.5 * ( h(:,j,k) + h(:,j,k-1) ) * GV%H_to_m ! Depth of center of layer K
        ! Mixed-layer depth, using sigma-0 (surface reference pressure)
        deltaRhoAtKm1(:) = deltaRhoAtK(:) ! Store value from previous iteration of K
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pRef_MLD, deltaRhoAtK, is-1, ie-is+3, tv%eqn_of_state)
        deltaRhoAtK(:) = deltaRhoAtK(:) - rhoSurf(:) ! Density difference between layer K and surface
        do i = is-1, ie+1
          ddRho = deltaRhoAtK(i) - deltaRhoAtKm1(i)
          if ((CS%MLD(i,j)==0.) .and. (ddRho>0.) .and. &
              (deltaRhoAtKm1(i)<CS%MLE_density_diff) .and. (deltaRhoAtK(i)>=CS%MLE_density_diff)) then
            aFac = ( CS%MLE_density_diff - deltaRhoAtKm1(i) ) / ddRho
            CS%MLD(i,j) = dK(i) * aFac + dKm1(i) * (1. - aFac)
          endif
        enddo ! i-loop
      enddo ! k-loop
      do i = is-1, ie+1
        if ((CS%MLD(i,j)==0.) .and. (deltaRhoAtK(i)<CS%MLE_density_diff)) CS%MLD(i,j) = dK(i) ! Assume mixing to the bottom
      enddo
    enddo ! j-loop
  elseif (CS%MLE_use_PBL_MLD .and. CS%MLE_MLD_decay_time>=0.) then
    if (.not. associated(MLD)) call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "Argument MLD was not associated!")
    if (CS%debug) then
      call hchksum(CS%MLD_filtered,'mixed_layer_restrat: MLD_filtered',G%HI,haloshift=1)
      call hchksum(MLD,'mixed_layer_restrat: MLD in',G%HI,haloshift=1)
    endif
    aFac = CS%MLE_MLD_decay_time / ( dt + CS%MLE_MLD_decay_time )
    bFac = dt / ( dt + CS%MLE_MLD_decay_time )
    do j = js-1, je+1 ; do i = is-1, ie+1
      ! Expression bFac*MLD(i,j) + aFac*CS%MLD_filtered(i,j) is the time-filtered
      ! (running mean) of MLD. The max() allows the "running mean" to be reset
      ! instantly to a deeper MLD.
      CS%MLD_filtered(i,j) = max( MLD(i,j), bFac*MLD(i,j) + aFac*CS%MLD_filtered(i,j) )
      CS%MLD(i,j) = CS%MLD_filtered(i,j)
    enddo ; enddo
  else
    call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "No MLD to use for MLE parameterization.")
  endif

  uDml(:) = 0.0 ; vDml(:) = 0.0
  I4dt = 0.25 / dt
  g_Rho0 = GV%g_Earth/GV%Rho0
  h_neglect = GV%H_subroundoff
  dz_neglect = GV%H_subroundoff*GV%H_to_m

  p0(:) = 0.0
!$OMP parallel default(none) shared(is,ie,js,je,G,GV,htot,Rml_av,tv,p0,h,h_avail,         &
!$OMP                               h_neglect,g_Rho0,I4dt,CS,uhml,uhtr,dt,vhml,vhtr,   &
!$OMP                               utimescale_diag,vtimescale_diag,fluxes,dz_neglect, &
!$OMP                               nz,MLD,uDml_diag,vDml_diag)                        &
!$OMP                       private(Rho0,h_vel,u_star,absf,mom_mixrate,timescale,      &
!$OMP                               a,IhTot,zIHbelowVel,hAtVel,zIHaboveVel)            &
!$OMP                       firstprivate(uDml,vDml)
!$OMP do
  do j=js-1,je+1
    do i=is-1,ie+1
      htot(i,j) = 0.0 ; Rml_av(i,j) = 0.0
    enddo
    do k=1,nz
      call calculate_density(tv%T(:,j,k),tv%S(:,j,k),p0,Rho0(:),is-1,ie-is+3,tv%eqn_of_state)
      do i=is-1,ie+1
        if (htot(i,j) < CS%MLD(i,j)) then
          Rml_av(i,j) = Rml_av(i,j) + h(i,j,k)*Rho0(i)
          htot(i,j) = htot(i,j) + h(i,j,k)
        endif
        h_avail(i,j,k) = max(I4dt*G%areaT(i,j)*(h(i,j,k)-GV%Angstrom),0.0)
      enddo
    enddo

    do i=is-1,ie+1
      Rml_av(i,j) = -(g_Rho0*Rml_av(i,j)) / (htot(i,j) + h_neglect)
    enddo
  enddo

  if (CS%debug) then
    call hchksum(h,'mixed_layer_restrat: h',G%HI,haloshift=1)
    call hchksum(fluxes%ustar,'mixed_layer_restrat: u*',G%HI,haloshift=1)
    call hchksum(CS%MLD,'mixed_layer_restrat: MLD',G%HI,haloshift=1)
    call hchksum(Rml_av,'mixed_layer_restrat: rml',G%HI,haloshift=1)
  endif

! TO DO:
!   1. Mixing extends below the mixing layer to the mixed layer.  Find it!
!   2. Add exponential tail to streamfunction?

!   U - Component
!$OMP do
  do j=js,je ; do I=is-1,ie
    h_vel = 0.5*((htot(i,j) + htot(i+1,j)) + h_neglect) * GV%H_to_m

    u_star = 0.5*(fluxes%ustar(i,j) + fluxes%ustar(i+1,j))
    absf = 0.5*(abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I,J)))
    ! peak ML visc: u_star * 0.41 * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
    ! momentum mixing rate: pi^2*visc/h_ml^2
    ! 0.41 is the von Karmen constant, 9.8696 = pi^2.
    mom_mixrate = (0.41*9.8696)*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+dz_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)

    timescale = timescale * CS%ml_restrat_coef
!     timescale = timescale*(2?)*(L_def/L_MLI)*min(EKE/MKE,1.0 + G%dyCv(i,j)**2/L_def**2))


    uDml(I) = timescale * G%mask2dCu(I,j)*G%dyCu(I,j)* &
        G%IdxCu(I,j)*(Rml_av(i+1,j)-Rml_av(i,j)) * (h_vel**2 * GV%m_to_H)

    if (uDml(i) == 0) then
      do k=1,nz ; uhml(I,j,k) = 0.0 ; enddo
    else
      IhTot = 2.0 / ((htot(i,j) + htot(i+1,j)) + h_neglect)
      zIHbelowVel = 0.0
      ! a(k) relates the sublayer transport to uDml with a linear profile.
      ! The sum of a(k) through the mixed layers must be 0.
      do k=1,nz
        hAtVel = 0.5*(h(i,j,k) + h(i+1,j,k))
        zIHaboveVel = zIHbelowVel                    ! z/H for upper interface
        zIHbelowVel = zIHbelowVel - (hAtVel * IhTot) ! z/H for lower interface
        a(k) = PSI( zIHaboveVel ) - PSI( zIHbelowVel )
        if (a(k)*uDml(I) > 0.0) then
          if (a(k)*uDml(I) > h_avail(i,j,k)) uDml(I) = h_avail(i,j,k) / a(k)
        elseif (a(k)*uDml(I) < 0.0) then
          if (-a(k)*uDml(I) > h_avail(i+1,j,k)) uDml(I) = -h_avail(i+1,j,k) / a(k)
        endif
      enddo
      do k=1,nz
        uhml(I,j,k) = a(k)*uDml(I)
        uhtr(I,j,k) = uhtr(I,j,k) + uhml(I,j,k)*dt
      enddo
    endif

    utimescale_diag(I,j) = timescale
    uDml_diag(I,j) = uDml(I)
  enddo ; enddo

!  V- component
!$OMP do
  do J=js-1,je ; do i=is,ie
    h_vel = 0.5*((htot(i,j) + htot(i,j+1)) + h_neglect) * GV%H_to_m

    u_star = 0.5*(fluxes%ustar(i,j) + fluxes%ustar(i,j+1))
    absf = 0.5*(abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))
    ! peak ML visc: u_star * 0.41 * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
    ! momentum mixing rate: pi^2*visc/h_ml^2
    ! 0.41 is the von Karmen constant, 9.8696 = pi^2.
    mom_mixrate = (0.41*9.8696)*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+dz_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)

    timescale = timescale * CS%ml_restrat_coef
!     timescale = timescale*(2?)*(L_def/L_MLI)*min(EKE/MKE,1.0 + G%dyCv(i,j)**2/L_def**2))

    vDml(i) = timescale * G%mask2dCv(i,J)*G%dxCv(i,J)* &
        G%IdyCv(i,J)*(Rml_av(i,j+1)-Rml_av(i,j)) * (h_vel**2 * GV%m_to_H)
    if (vDml(i) == 0) then
      do k=1,nz ; vhml(i,J,k) = 0.0 ; enddo
    else
      IhTot = 2.0 / ((htot(i,j) + htot(i,j+1)) + h_neglect)
      zIHbelowVel = 0.0
      ! a(k) relates the sublayer transport to uDml with a linear profile.
      ! The sum of a(k) through the mixed layers must be 0.
      do k=1,nz
        hAtVel = 0.5*(h(i,j,k) + h(i,j+1,k))
        zIHaboveVel = zIHbelowVel                    ! z/H for upper interface
        zIHbelowVel = zIHbelowVel - (hAtVel * IhTot) ! z/H for lower interface
        a(k) = PSI( zIHaboveVel ) - PSI( zIHbelowVel )
        if (a(k)*vDml(i) > 0.0) then
          if (a(k)*vDml(i) > h_avail(i,j,k)) vDml(i) = h_avail(i,j,k) / a(k)
        elseif (a(k)*vDml(i) < 0.0) then
          if (-a(k)*vDml(i) > h_avail(i,j+1,k)) vDml(i) = -h_avail(i,j+1,k) / a(k)
        endif
      enddo
      do k=1,nz
        vhml(i,J,k) = a(k)*vDml(i)
        vhtr(i,J,k) = vhtr(i,J,k) + vhml(i,J,k)*dt
      enddo
    endif

    vtimescale_diag(i,J) = timescale
    vDml_diag(i,J) = vDml(i)
  enddo ; enddo

!$OMP do
  do j=js,je ; do k=1,nz ; do i=is,ie
    h(i,j,k) = h(i,j,k) - dt*G%IareaT(i,j) * &
        ((uhml(I,j,k) - uhml(I-1,j,k)) + (vhml(i,J,k) - vhml(i,J-1,k)))
  enddo ; enddo ; enddo
!$OMP end parallel

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  ! This needs to happen after the H update and before the next post_data.
  call diag_update_remap_grids(CS%diag)

  ! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_urestrat_time > 0) call post_data(CS%id_urestrat_time, utimescale_diag, CS%diag)
    if (CS%id_vrestrat_time > 0) call post_data(CS%id_vrestrat_time, vtimescale_diag, CS%diag)
    if (CS%id_uhml          > 0) call post_data(CS%id_uhml, uhml, CS%diag)
    if (CS%id_vhml          > 0) call post_data(CS%id_vhml, vhml, CS%diag)
    if (CS%id_MLD           > 0) call post_data(CS%id_MLD, CS%MLD, CS%diag)
    if (CS%id_Rml           > 0) call post_data(CS%id_Rml, Rml_av, CS%diag)
    if (CS%id_uDml          > 0) call post_data(CS%id_uDml, uDml_diag, CS%diag)
    if (CS%id_vDml          > 0) call post_data(CS%id_vDml, vDml_diag, CS%diag)

    if (CS%id_uml > 0) then
      do J=js,je ; do i=is,ie
        h_vel = 0.5*((htot(i,j) + htot(i+1,j)) + h_neglect)
        uDml_diag(I,j) = uDml_diag(I,j) / (0.01*h_vel) * G%IdyCu(I,j) * (PSI(0.)-PSI(-.01))
      enddo ; enddo
      call post_data(CS%id_uml, uDml_diag, CS%diag)
    endif
    if (CS%id_vml > 0) then
      do J=js,je ; do i=is,ie
        h_vel = 0.5*((htot(i,j) + htot(i,j+1)) + h_neglect)
        vDml_diag(i,J) = vDml_diag(i,J) / (0.01*h_vel) * G%IdxCv(i,J) * (PSI(0.)-PSI(-.01))
      enddo ; enddo
      call post_data(CS%id_vml, vDml_diag, CS%diag)
    endif
  endif

end subroutine mixedlayer_restrat_general


!> This subroutine does interface depth diffusion.  The fluxes are
!! limited to give positive definiteness, and the diffusivities are
!! limited to guarantee stability.
subroutine mixedlayer_restrat_BML(h, uhtr, vhtr, tv, fluxes, dt, G, GV, CS)
  type(ocean_grid_type),                     intent(in)    :: G      !< ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h      !< layer thickness (H units = m or kg/m2)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uhtr   !< accumulated zonal mass flux (m3 or kg)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vhtr   !< accumulated merid mass flux (m3 or kg)
  type(thermo_var_ptrs),                     intent(in)    :: tv     !< thermodynamic variables structure
  type(forcing),                             intent(in)    :: fluxes !< pointers to forcing fields
  real,                                      intent(in)    :: dt     !< time increment (sec)
  type(mixedlayer_restrat_CS),               pointer       :: CS     !< module control structure

  real :: uhml(SZIB_(G),SZJ_(G),SZK_(G)) ! zonal mixed layer transport (m3/s or kg/s)
  real :: vhml(SZI_(G),SZJB_(G),SZK_(G)) ! merid mixed layer transport (m3/s or kg/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    h_avail               ! The volume available for diffusion out of each face of each
                          ! sublayer of the mixed layer, divided by dt, in units
                          ! of H m2 s-1 (i.e., m3 s-1 or kg s-1).
  real, dimension(SZI_(G),SZJ_(G)) :: &
    htot, &               ! The sum of the thicknesses of layers in the mixed layer (H units)
    Rml_av                ! g_Rho0 times the average mixed layer density (m s-2)
  real :: g_Rho0          ! G_Earth/Rho0 (m4 s-2 kg-1)
  real :: Rho0(SZI_(G))   ! Potential density relative to the surface (kg m-3)
  real :: p0(SZI_(G))     ! A pressure of 0 (Pa)

  real :: h_vel           ! htot interpolated onto velocity points (meter; not H)
  real :: absf            ! absolute value of f, interpolated to velocity points (s-1)
  real :: u_star          ! surface friction velocity, interpolated to velocity points (m s-1)
  real :: mom_mixrate     ! rate at which momentum is homogenized within mixed layer (s-1)
  real :: timescale       ! mixing growth timescale (sec)
  real :: h_neglect       ! tiny thickness usually lost in roundoff and can be neglected (H units)
  real :: dz_neglect      ! tiny thickness (in m) that usually lost in roundoff and can be neglected (meter)
  real :: I4dt            ! 1/(4 dt)
  real :: I2htot          ! Twice the total mixed layer thickness at velocity points (H units)
  real :: z_topx2         ! depth of the top of a layer at velocity points (H units)
  real :: hx2             ! layer thickness at velocity points (H units)
  real :: a(SZK_(G))      ! A nondimensional value relating the overall flux
                          ! magnitudes (uDml & vDml) to the realized flux in a
                          ! layer.  The vertical sum of a() through the pieces of
                          ! the mixed layer must be 0.
  real :: uDml(SZIB_(G))  ! The zonal and meridional volume fluxes in the upper
  real :: vDml(SZI_(G))   ! half of the mixed layer in H m2 s-1 (m3 s-1 or kg s-1).
  real :: utimescale_diag(SZIB_(G),SZJ_(G)) ! The restratification timescales
  real :: vtimescale_diag(SZI_(G),SZJB_(G)) ! in the zonal and meridional
                                            ! directions (sec), stored in 2-D
                                            ! arrays for diagnostic purposes.
  real :: uDml_diag(SZIB_(G),SZJ_(G)), vDml_diag(SZI_(G),SZJB_(G))
  logical :: use_EOS    ! If true, density is calculated from T & S using an equation of state.

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkml
  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nkml = GV%nkml

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "Module must be initialized before it is used.")
  if ((nkml<2) .or. (CS%ml_restrat_coef<=0.0)) return

  uDml(:)    = 0.0 ; vDml(:) = 0.0
  I4dt       = 0.25 / dt
  g_Rho0     = GV%g_Earth/GV%Rho0
  use_EOS    = associated(tv%eqn_of_state)
  h_neglect  = GV%H_subroundoff
  dz_neglect = GV%H_subroundoff*GV%H_to_m

  if (.not.use_EOS) call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "An equation of state must be used with this module.")

  ! Fix this later for nkml >= 3.

  p0(:) = 0.0
!$OMP parallel default(none) shared(is,ie,js,je,G,GV,htot,Rml_av,tv,p0,h,h_avail,      &
!$OMP                               h_neglect,g_Rho0,I4dt,CS,uhml,uhtr,dt,vhml,vhtr,   &
!$OMP                               utimescale_diag,vtimescale_diag,fluxes,dz_neglect, &
!$OMP                               uDml_diag,vDml_diag,nkml)                          &
!$OMP                       private(Rho0,h_vel,u_star,absf,mom_mixrate,timescale,      &
!$OMP                               I2htot,z_topx2,hx2,a)                              &
!$OMP                       firstprivate(uDml,vDml)
!$OMP do
  do j=js-1,je+1
    do i=is-1,ie+1
      htot(i,j) = 0.0 ; Rml_av(i,j) = 0.0
    enddo
    do k=1,nkml
      call calculate_density(tv%T(:,j,k),tv%S(:,j,k),p0,Rho0(:),is-1,ie-is+3,tv%eqn_of_state)
      do i=is-1,ie+1
        Rml_av(i,j) = Rml_av(i,j) + h(i,j,k)*Rho0(i)
        htot(i,j) = htot(i,j) + h(i,j,k)
        h_avail(i,j,k) = max(I4dt*G%areaT(i,j)*(h(i,j,k)-GV%Angstrom),0.0)
      enddo
    enddo

    do i=is-1,ie+1
      Rml_av(i,j) = (g_Rho0*Rml_av(i,j)) / (htot(i,j) + h_neglect)
    enddo
  enddo

! TO DO:
!   1. Mixing extends below the mixing layer to the mixed layer.  Find it!
!   2. Add exponential tail to streamfunction?

!   U - Component
!$OMP do
  do j=js,je
    do i=is,ie ; utimescale_diag(i,j) = 0.0 ; enddo
    do i=is,ie ; vtimescale_diag(i,j) = 0.0 ; enddo
    do I=is-1,ie
      h_vel = 0.5*(htot(i,j) + htot(i+1,j)) * GV%H_to_m

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
          G%IdxCu(I,j)*(Rml_av(i+1,j)-Rml_av(i,j)) * (h_vel**2 * GV%m_to_H)

      if (uDml(i) == 0) then
        do k=1,nkml ; uhml(I,j,k) = 0.0 ; enddo
      else
        I2htot = 1.0 / (htot(i,j) + htot(i+1,j) + h_neglect)
        z_topx2 = 0.0
        ! a(k) relates the sublayer transport to uDml with a linear profile.
        ! The sum of a(k) through the mixed layers must be 0.
        do k=1,nkml
          hx2 = (h(i,j,k) + h(i+1,j,k) + h_neglect)
          a(k) = (hx2 * I2htot) * (2.0 - 4.0*(z_topx2+0.5*hx2)*I2htot)
          z_topx2 = z_topx2 + hx2
          if (a(k)*uDml(I) > 0.0) then
            if (a(k)*uDml(I) > h_avail(i,j,k)) uDml(I) = h_avail(i,j,k) / a(k)
          else
            if (-a(k)*uDml(I) > h_avail(i+1,j,k)) uDml(I) = -h_avail(i+1,j,k)/a(k)
          endif
        enddo
        do k=1,nkml
          uhml(I,j,k) = a(k)*uDml(I)
          uhtr(I,j,k) = uhtr(I,j,k) + uhml(I,j,k)*dt
        enddo
      endif
    enddo
    uDml_diag(is:ie,j) = uDml(is:ie)
  enddo

!  V- component
!$OMP do
  do J=js-1,je ; do i=is,ie
    h_vel = 0.5*(htot(i,j) + htot(i,j+1)) * GV%H_to_m

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
        G%IdyCv(i,J)*(Rml_av(i,j+1)-Rml_av(i,j)) * (h_vel**2 * GV%m_to_H)
    if (vDml(i) == 0) then
      do k=1,nkml ; vhml(i,J,k) = 0.0 ; enddo
    else
      I2htot = 1.0 / (htot(i,j) + htot(i,j+1) + h_neglect)
      z_topx2 = 0.0
      ! a(k) relates the sublayer transport to uDml with a linear profile.
      ! The sum of a(k) through the mixed layers must be 0.
      do k=1,nkml
        hx2 = (h(i,j,k) + h(i,j+1,k) + h_neglect)
        a(k) = (hx2 * I2htot) * (2.0 - 4.0*(z_topx2+0.5*hx2)*I2htot)
        z_topx2 = z_topx2 + hx2
        if (a(k)*vDml(i) > 0.0) then
          if (a(k)*vDml(i) > h_avail(i,j,k)) vDml(i) = h_avail(i,j,k) / a(k)
        else
          if (-a(k)*vDml(i) > h_avail(i,j+1,k)) vDml(i) = -h_avail(i,j+1,k)/a(k)
        endif
      enddo
      do k=1,nkml
        vhml(i,J,k) = a(k)*vDml(i)
        vhtr(i,J,k) = vhtr(i,J,k) + vhml(i,J,k)*dt
      enddo
    endif
    enddo
    vDml_diag(is:ie,j) = vDml(is:ie)
  enddo

!$OMP do
  do j=js,je ; do k=1,nkml ; do i=is,ie
    h(i,j,k) = h(i,j,k) - dt*G%IareaT(i,j) * &
        ((uhml(I,j,k) - uhml(I-1,j,k)) + (vhml(i,J,k) - vhml(i,J-1,k)))
  enddo ; enddo ; enddo
!$OMP end parallel

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  call diag_update_remap_grids(CS%diag)

  ! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag) .and. &
    ((CS%id_urestrat_time > 0)  .or. (CS%id_vrestrat_time > 0))) then
    call post_data(CS%id_urestrat_time, utimescale_diag, CS%diag)
    call post_data(CS%id_vrestrat_time, vtimescale_diag, CS%diag)
  endif
  if (query_averaging_enabled(CS%diag) .and. &
      ((CS%id_uhml>0) .or. (CS%id_vhml>0))) then
    do k=nkml+1,nz
      do j=js,je ; do I=Isq,Ieq ; uhml(I,j,k) = 0.0 ; enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie ; vhml(i,J,k) = 0.0 ; enddo ; enddo
    enddo
    if (CS%id_uhml > 0) call post_data(CS%id_uhml, uhml,      CS%diag)
    if (CS%id_vhml > 0) call post_data(CS%id_vhml, vhml,      CS%diag)
    if (CS%id_MLD  > 0) call post_data(CS%id_MLD,  htot,      CS%diag)
    if (CS%id_Rml  > 0) call post_data(CS%id_Rml,  Rml_av,    CS%diag)
    if (CS%id_uDml > 0) call post_data(CS%id_uDml, uDml_diag, CS%diag)
    if (CS%id_vDml > 0) call post_data(CS%id_vDml, vDml_diag, CS%diag)
  endif

end subroutine mixedlayer_restrat_BML


!> Initialize the mixedlayer restratification module
logical function mixedlayer_restrat_init(Time, G, GV, param_file, diag, CS)
  type(time_type),             intent(in)    :: Time       !< current model time
  type(ocean_grid_type),       intent(inout) :: G          !< ocean grid structure
  type(verticalGrid_type),     intent(in)    :: GV         !< ocean vertical grid structure
  type(param_file_type),       intent(in)    :: param_file !< parameter file to parse
  type(diag_ctrl), target,     intent(inout) :: diag       !< regulate diagnostics
  type(mixedlayer_restrat_CS), pointer       :: CS         !< module control structure

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=48)  :: flux_units

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "MIXEDLAYER_RESTRAT", mixedlayer_restrat_init, &
             "If true, a density-gradient dependent re-stratifying \n"//&
             "flow is imposed in the mixed layer. Can be used in ALE mode\n"//&
             "without restriction but in layer mode can only be used if\n"//&
             "BULKMIXEDLAYER is true.", default=.false.)
  if (.not. mixedlayer_restrat_init) return

  if (.not.associated(CS)) then
    call MOM_error(FATAL, "mixedlayer_restrat_init called without an "// &
                            "associated control structure.")
  endif

  ! Nonsense values to cause problems when these parameters are not used
  CS%MLE_MLD_decay_time = -9.e9
  CS%MLE_density_diff = -9.e9
  CS%MLE_tail_dh = -9.e9
  CS%MLE_use_PBL_MLD = .false.

  call get_param(param_file, mod, "DEBUG", CS%debug, default=.false., do_not_log=.true.)
  call get_param(param_file, mod, "FOX_KEMPER_ML_RESTRAT_COEF", CS%ml_restrat_coef, &
             "A nondimensional coefficient that is proportional to \n"//&
             "the ratio of the deformation radius to the dominant \n"//&
             "lengthscale of the submesoscale mixed layer \n"//&
             "instabilities, times the minimum of the ratio of the \n"//&
             "mesoscale eddy kinetic energy to the large-scale \n"//&
             "geostrophic kinetic energy or 1 plus the square of the \n"//&
             "grid spacing over the deformation radius, as detailed \n"//&
             "by Fox-Kemper et al. (2010)", units="nondim", default=0.0)
  ! We use GV%nkml to distinguish between the old and new implementation of MLE.
  ! The old implementation only works for the layer model with nkml>0.
  if (GV%nkml==0) then
    call get_param(param_file, mod, "MLE_USE_PBL_MLD", CS%MLE_use_PBL_MLD, &
             "If true, the MLE parameterization will use the mixed-layer\n"//&
             "depth provided by the active PBL parameterization. If false,\n"//&
             "MLE will estimate a MLD based on a density difference with the\n"//&
             "surface using the parameter MLE_DENSITY_DIFF.", default=.false.)
    if (CS%MLE_use_PBL_MLD) then
      call get_param(param_file, mod, "MLE_MLD_DECAY_TIME", CS%MLE_MLD_decay_time, &
             "When MLE_USE_PBL_MLD is true, the PBL provided active miing layer\n"//&
             "depth is running mean average with this time-scale when the MLD\n"//&
             "is retreating. When the MLD deepens below the current running-mean\n"//&
             "the runn-mean is instantaneously set to the current MLD.", units="s", default=0.)
    else
      call get_param(param_file, mod, "MLE_DENSITY_DIFF", CS%MLE_density_diff, &
             "Density difference used to detect the mixed-layer\n"//&
             "depth used for the mixed-layer eddy parameterization\n"//&
             "by Fox-Kemper et al. (2010)", units="kg/m3", default=0.03)
    endif
    call get_param(param_file, mod, "MLE_TAIL_DH", CS%MLE_tail_dh, &
             "Fraction by which to extend the mixed-layer restratification\n"//&
             "depth used for a smoother stream function at the base of\n"//&
             "the mixed-layer.", units="nondim", default=0.0)
  endif

  CS%diag => diag

  if (GV%Boussinesq) then ; flux_units = "meter3 second-1"
  else ; flux_units = "kilogram second-1" ; endif

  CS%id_uhml = register_diag_field('ocean_model', 'uhml', diag%axesCuL, Time, &
      'Zonal Thickness Flux to Restratify Mixed Layer', flux_units, &
      y_cell_method='sum', v_extensive=.true.)
  CS%id_vhml = register_diag_field('ocean_model', 'vhml', diag%axesCvL, Time, &
      'Meridional Thickness Flux to Restratify Mixed Layer', flux_units, &
      x_cell_method='sum', v_extensive=.true.)
  CS%id_urestrat_time = register_diag_field('ocean_model', 'MLu_restrat_time', diag%axesCu1, Time, &
      'Mixed Layer Zonal Restratification Timescale', 'second')
  CS%id_vrestrat_time = register_diag_field('ocean_model', 'MLv_restrat_time', diag%axesCv1, Time, &
      'Mixed Layer Meridional Restratification Timescale', 'second')
  CS%id_MLD = register_diag_field('ocean_model', 'MLD_restrat', diag%axesT1, Time, &
      'Mixed Layer Depth as used in the mixed-layer restratification parameterization', 'meter')
  CS%id_Rml = register_diag_field('ocean_model', 'ML_buoy_restrat', diag%axesT1, Time, &
      'Mixed Layer Buoyancy as used in the mixed-layer restratification parameterization', 'm/s^2')
  CS%id_uDml = register_diag_field('ocean_model', 'udml_restrat', diag%axesCu1, Time, &
      'Transport stream function amplitude for zonal restratification of mixed layer', 'm3/s')
  CS%id_vDml = register_diag_field('ocean_model', 'vdml_restrat', diag%axesCv1, Time, &
      'Transport stream function amplitude for meridional restratification of mixed layer', 'm3/s')
  CS%id_uml = register_diag_field('ocean_model', 'uml_restrat', diag%axesCu1, Time, &
      'Surface zonal velocity component of mixed layer restratification', 'm/s')
  CS%id_vml = register_diag_field('ocean_model', 'vml_restrat', diag%axesCv1, Time, &
      'Surface meridional velocity component of mixed layer restratification', 'm/s')

  ! If MLD_filtered is being used, we need to update halo regions after a restart
  if (associated(CS%MLD_filtered)) call pass_var(CS%MLD_filtered, G%domain)

end function mixedlayer_restrat_init

!> Allocate and regsiter fields in the mixedlayer restratification structure for restarts
subroutine mixedlayer_restrat_register_restarts(HI, param_file, CS, restart_CS)
  type(hor_index_type),        intent(in)    :: HI         !< Horizontal index structure
  type(param_file_type),       intent(in)    :: param_file !< Parameter file to parse
  type(mixedlayer_restrat_CS), pointer       :: CS         !< Module control structure
  type(MOM_restart_CS),        pointer       :: restart_CS !< Restart structure
  ! Local variables
  type(vardesc) :: vd
  logical :: mixedlayer_restrat_init

  ! Check to see if this module will be used
  call get_param(param_file, mod, "MIXEDLAYER_RESTRAT", mixedlayer_restrat_init, &
             default=.false., do_not_log=.true.)
  if (.not. mixedlayer_restrat_init) return

  ! Allocate the control structure. CS will be later populated by mixedlayer_restrat_init()
  if (associated(CS)) call MOM_error(FATAL, &
       "mixedlayer_restrat_register_restarts called with an associated control structure.")
  allocate(CS)

  ! CS%MLD is used either for the internally diagnosed MLD or
  ! for keep a running mean of the PBL's actively mixed MLD.
  allocate(CS%MLD(HI%isd:HI%ied,HI%jsd:HI%jed)) ; CS%MLD(:,:) = 0.

  call get_param(param_file, mod, "MLE_USE_PBL_MLD", CS%MLE_use_PBL_MLD, &
             default=.false., do_not_log=.true.)
  if (CS%MLE_use_PBL_MLD) then
    allocate(CS%MLD_filtered(HI%isd:HI%ied,HI%jsd:HI%jed)) ; CS%MLD_filtered(:,:) = 0.
    vd = var_desc("MLD_MLE_filtered","m","Time-filtered MLD for use in MLE", &
                  hor_grid='h', z_grid='1')
    call register_restart_field(CS%MLD_filtered, vd, .false., restart_CS)

  endif

end subroutine mixedlayer_restrat_register_restarts

!> \namespace mom_mixed_layer_restrat
!!
!!    The subroutine in this file implements a parameterization of
!!  unresolved viscous mixed layer restratification of the mixed layer
!!  as described in Fox-Kemper, Ferrari and Hallberg (JPO, 2008), and
!!  whose impacts are described in Fox-Kemper et al. (Ocean Modelling,
!!  2011).  This is derived in part from the older parameterizaton
!!  that is described in Hallberg (Aha Hulikoa, 2003), which this new
!!  parameterization surpasses, which in turn is based on the
!!  subinertial mixed layer theory of Young (JPO, 1994).  There is no
!!  net horizontal volume transport due to this parameterization, and
!!  no direct effect below the mixed layer.
!!
!!  This parameterization sets the restratification timescale to agree
!!  high-resolution studies of mixed layer restratification.  The run-time
!!  parameter FOX_KEMPER_ML_RESTRAT_COEF is a nondimensional number
!!  of order a few tens, proportional to the ratio of the deformation
!!  radius or the gridscale (whichever is smaller to the dominant
!!  horizontal lengthscale of the submesoscale mixed layer
!!  instabilities.
!!
!!  Macros written all in capital letters are defined in MOM_memory.h.
!!
!!  \section section_gridlayout MOM grid layout
!!
!!  A small fragment of the grid is shown below:
!!
!! \verbatim
!!    j+1  x ^ x ^ x
!!
!!    j+1  > o > o >
!!
!!    j    x ^ x ^ x
!!
!!    j    > o > o >
!!
!!    j-1  x ^ x ^ x
!!
!!        i-1  i  i+1
!!
!!           i  i+1
!!
!! \endverbatim
!!
!!  Fields at each point
!!  * x =  q, CoriolisBu
!!  * ^ =  v, PFv, CAv, vh, diffv, tauy, vbt, vhtr
!!  * > =  u, PFu, CAu, uh, diffu, taux, ubt, uhtr
!!  * o =  h, bathyT, eta, T, S, tr
!!
!!  The boundaries always run through q grid points (x).
!!


end module MOM_mixed_layer_restrat
