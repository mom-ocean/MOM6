!> \brief Parameterization of mixed layer restratification by unresolved mixed-layer eddies.
module MOM_mixed_layer_restrat

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,     only : hchksum
use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_diag_mediator, only : diag_update_remap_grids
use MOM_domains,       only : pass_var, To_West, To_South, Omit_Corners
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_forcing_type,  only : mech_forcing
use MOM_grid,          only : ocean_grid_type
use MOM_hor_index,     only : hor_index_type
use MOM_io,            only : vardesc, var_desc
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_restart,       only : register_restart_field, query_initialized, MOM_restart_CS
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_EOS,           only : calculate_density

implicit none ; private

#include <MOM_memory.h>

public mixedlayer_restrat
public mixedlayer_restrat_init
public mixedlayer_restrat_register_restarts

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for mom_mixed_layer_restrat
type, public :: mixedlayer_restrat_CS ; private
  real    :: ml_restrat_coef       !< A non-dimensional factor by which the instability is enhanced
                                   !! over what would be predicted based on the resolved gradients
                                   !! [nondim].  This increases with grid spacing^2, up to something
                                   !! of order 500.
  real    :: ml_restrat_coef2      !< As for ml_restrat_coef but using the slow filtered MLD [nondim].
  real    :: front_length          !< If non-zero, is the frontal-length scale [L ~> m] used to calculate the
                                   !! upscaling of buoyancy gradients that is otherwise represented
                                   !! by the parameter FOX_KEMPER_ML_RESTRAT_COEF. If MLE_FRONT_LENGTH is
                                   !! non-zero, it is recommended to set FOX_KEMPER_ML_RESTRAT_COEF=1.0.
  logical :: MLE_use_PBL_MLD       !< If true, use the MLD provided by the PBL parameterization.
                                   !! if false, MLE will calculate a MLD based on a density difference
                                   !! based on the parameter MLE_DENSITY_DIFF.
  real    :: MLE_MLD_decay_time    !< Time-scale to use in a running-mean when MLD is retreating [T ~> s].
  real    :: MLE_MLD_decay_time2   !< Time-scale to use in a running-mean when filtered MLD is retreating [T ~> s].
  real    :: MLE_density_diff      !< Density difference used in detecting mixed-layer depth [R ~> kg m-3].
  real    :: MLE_tail_dh           !< Fraction by which to extend the mixed-layer restratification
                                   !! depth used for a smoother stream function at the base of
                                   !! the mixed-layer [nondim].
  real    :: MLE_MLD_stretch       !< A scaling coefficient for stretching/shrinking the MLD used in
                                   !! the MLE scheme [nondim]. This simply multiplies MLD wherever used.
  logical :: MLE_use_MLD_ave_bug   !< If true, do not account for MLD mismatch to interface positions.
  logical :: debug = .false.       !< If true, calculate checksums of fields for debugging.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  real, dimension(:,:), pointer :: &
         MLD_filtered => NULL(), &   !< Time-filtered MLD [H ~> m or kg m-2]
         MLD_filtered_slow => NULL() !< Slower time-filtered MLD [H ~> m or kg m-2]

  !>@{
  !! Diagnostic identifier
  integer :: id_urestrat_time = -1
  integer :: id_vrestrat_time = -1
  integer :: id_uhml = -1
  integer :: id_vhml = -1
  integer :: id_MLD  = -1
  integer :: id_Rml  = -1
  integer :: id_uDml = -1
  integer :: id_vDml = -1
  integer :: id_uml  = -1
  integer :: id_vml  = -1
  !>@}

end type mixedlayer_restrat_CS

character(len=40)  :: mdl = "MOM_mixed_layer_restrat" !< This module's name.

contains

!> Driver for the mixed-layer restratification parameterization.
!! The code branches between two different implementations depending
!! on whether the bulk-mixed layer or a general coordinate are in use.
subroutine mixedlayer_restrat(h, uhtr, vhtr, tv, forces, dt, MLD, VarMix, G, GV, US, CS)
  type(ocean_grid_type),                     intent(inout) :: G      !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                     !! [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                     !! [H L2 ~> m3 or kg]
  type(thermo_var_ptrs),                     intent(in)    :: tv     !< Thermodynamic variables structure
  type(mech_forcing),                        intent(in)    :: forces !< A structure with the driving mechanical forces
  real,                                      intent(in)    :: dt     !< Time increment [T ~> s]
  real, dimension(:,:),                      pointer       :: MLD    !< Mixed layer depth provided by the
                                                                     !! PBL scheme [H ~> m or kg m-2]
  type(VarMix_CS),                           pointer       :: VarMix !< Container for derived fields
  type(mixedlayer_restrat_CS),               pointer       :: CS     !< Module control structure

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "Module must be initialized before it is used.")

  if (GV%nkml>0) then
    call mixedlayer_restrat_BML(h, uhtr, vhtr, tv, forces, dt, G, GV, US, CS)
  else
    call mixedlayer_restrat_general(h, uhtr, vhtr, tv, forces, dt, MLD, VarMix, G, GV, US, CS)
  endif

end subroutine mixedlayer_restrat

!> Calculates a restratifying flow in the mixed layer.
subroutine mixedlayer_restrat_general(h, uhtr, vhtr, tv, forces, dt, MLD_in, VarMix, G, GV, US, CS)
  ! Arguments
  type(ocean_grid_type),                     intent(inout) :: G      !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                     !!   [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                     !!   [H L2 ~> m3 or kg]
  type(thermo_var_ptrs),                     intent(in)    :: tv     !< Thermodynamic variables structure
  type(mech_forcing),                        intent(in)    :: forces !< A structure with the driving mechanical forces
  real,                                      intent(in)    :: dt     !< Time increment [T ~> s]
  real, dimension(:,:),                      pointer       :: MLD_in !< Mixed layer depth provided by the
                                                                     !! PBL scheme [m] (not H)
  type(VarMix_CS),                           pointer       :: VarMix !< Container for derived fields
  type(mixedlayer_restrat_CS),               pointer       :: CS     !< Module control structure
  ! Local variables
  real :: uhml(SZIB_(G),SZJ_(G),SZK_(G)) ! zonal mixed layer transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: vhml(SZI_(G),SZJB_(G),SZK_(G)) ! merid mixed layer transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    h_avail               ! The volume available for diffusion out of each face of each
                          ! sublayer of the mixed layer, divided by dt [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    MLD_fast, &           ! Mixed layer depth actually used in MLE restratification parameterization [H ~> m or kg m-2]
    htot_fast, &          ! The sum of the thicknesses of layers in the mixed layer [H ~> m or kg m-2]
    Rml_av_fast, &        ! g_Rho0 times the average mixed layer density [L2 Z-1 T-2 ~> m s-2]
    MLD_slow, &           ! Mixed layer depth actually used in MLE restratification parameterization [H ~> m or kg m-2]
    htot_slow, &          ! The sum of the thicknesses of layers in the mixed layer [H ~> m or kg m-2]
    Rml_av_slow           ! g_Rho0 times the average mixed layer density [L2 Z-1 T-2 ~> m s-2]
  real :: g_Rho0          ! G_Earth/Rho0 [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1]
  real :: rho_ml(SZI_(G)) ! Potential density relative to the surface [R ~> kg m-3]
  real :: p0(SZI_(G))     ! A pressure of 0 [Pa]

  real :: h_vel           ! htot interpolated onto velocity points [Z ~> m] (not H).
  real :: absf            ! absolute value of f, interpolated to velocity points [T-1 ~> s-1]
  real :: u_star          ! surface friction velocity, interpolated to velocity points [Z T-1 ~> m s-1].
  real :: mom_mixrate     ! rate at which momentum is homogenized within mixed layer [T-1 ~> s-1]
  real :: timescale       ! mixing growth timescale [T ~> s]
  real :: h_neglect       ! tiny thickness usually lost in roundoff so can be neglected [H ~> m or kg m-2]
  real :: dz_neglect      ! A tiny thickness that is usually lost in roundoff so can be neglected [Z ~> m]
  real :: I4dt            ! 1/(4 dt) [T-1 ~> s-1]
  real :: Ihtot,Ihtot_slow! Inverses of the total mixed layer thickness [H-1 ~> m-1 or m2 kg-1]
  real :: a(SZK_(G))      ! A non-dimensional value relating the overall flux
                          ! magnitudes (uDml & vDml) to the realized flux in a
                          ! layer.  The vertical sum of a() through the pieces of
                          ! the mixed layer must be 0.
  real :: b(SZK_(G))      ! As for a(k) but for the slow-filtered MLD
  real :: uDml(SZIB_(G))  ! The zonal and meridional volume fluxes in the upper
  real :: vDml(SZI_(G))   ! half of the mixed layer [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: uDml_slow(SZIB_(G))  ! The zonal and meridional volume fluxes in the upper
  real :: vDml_slow(SZI_(G))   ! half of the mixed layer [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: utimescale_diag(SZIB_(G),SZJ_(G)) ! restratification timescales in the zonal and
  real :: vtimescale_diag(SZI_(G),SZJB_(G)) ! meridional directions [T ~> s], stored in 2-D arrays
                                            ! for diagnostic purposes.
  real :: uDml_diag(SZIB_(G),SZJ_(G)), vDml_diag(SZI_(G),SZJB_(G))
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  real, dimension(SZI_(G)) :: rhoSurf, deltaRhoAtKm1, deltaRhoAtK ! Densities [R ~> kg m-3]
  real, dimension(SZI_(G)) :: dK, dKm1 ! Depths of layer centers [H ~> m or kg m-2].
  real, dimension(SZI_(G)) :: pRef_MLD ! A reference pressure for calculating the mixed layer densities [Pa].
  real, dimension(SZI_(G)) :: rhoAtK, rho1, d1, pRef_N2 ! Used for N2
  real :: aFac, bFac ! Nondimensional ratios [nondim]
  real :: ddRho    ! A density difference [R ~> kg m-3]
  real :: hAtVel, zpa, zpb, dh, res_scaling_fac
  real :: I_LFront ! The inverse of the frontal length scale [L-1 ~> m-1]
  logical :: proper_averaging, line_is_empty, keep_going, res_upscale

  real :: PSI, PSI1, z, BOTTOP, XP, DD ! For the following statement functions
  ! Stream function as a function of non-dimensional position within mixed-layer (F77 statement function)
  !PSI1(z) = max(0., (1. - (2.*z+1.)**2 ) )
  PSI1(z) = max(0., (1. - (2.*z+1.)**2 ) * (1. + (5./21.)*(2.*z+1.)**2) )
  BOTTOP(z) = 0.5*(1.-SIGN(1.,z+0.5)) ! =0 for z>-0.5, =1 for z<-0.5
  XP(z) = max(0., min(1., (-z-0.5)*2./(1.+2.*CS%MLE_tail_dh) ) )
  DD(z) = (1.-3.*(XP(z)**2)+2.*(XP(z)**3))**(1.+2.*CS%MLE_tail_dh)
  PSI(z) = max( PSI1(z), DD(z)*BOTTOP(z) ) ! Combines original PSI1 with tail

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not.associated(tv%eqn_of_state)) call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "An equation of state must be used with this module.")
  if (.not.associated(VarMix) .and. CS%front_length>0.) call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "The resolution argument, Rd/dx, was not associated.")

  if (CS%MLE_density_diff > 0.) then ! We need to calculate a mixed layer depth, MLD.
    !! TODO: use derivatives and mid-MLD pressure. Currently this is sigma-0. -AJA
    pRef_MLD(:) = 0.
    do j = js-1, je+1
      dK(:) = 0.5 * h(:,j,1) ! Depth of center of surface layer
      call calculate_density(tv%T(:,j,1), tv%S(:,j,1), pRef_MLD, rhoSurf, is-1, ie-is+3, &
                             tv%eqn_of_state, scale=US%kg_m3_to_R)
      deltaRhoAtK(:) = 0.
      MLD_fast(:,j) = 0.
      do k = 2, nz
        dKm1(:) = dK(:) ! Depth of center of layer K-1
        dK(:) = dK(:) + 0.5 * ( h(:,j,k) + h(:,j,k-1) ) ! Depth of center of layer K
        ! Mixed-layer depth, using sigma-0 (surface reference pressure)
        deltaRhoAtKm1(:) = deltaRhoAtK(:) ! Store value from previous iteration of K
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pRef_MLD, deltaRhoAtK, is-1, ie-is+3, &
                               tv%eqn_of_state, scale=US%kg_m3_to_R)
        do i = is-1,ie+1
          deltaRhoAtK(i) = deltaRhoAtK(i) - rhoSurf(i) ! Density difference between layer K and surface
        enddo
        do i = is-1, ie+1
          ddRho = deltaRhoAtK(i) - deltaRhoAtKm1(i)
          if ((MLD_fast(i,j)==0.) .and. (ddRho>0.) .and. &
              (deltaRhoAtKm1(i)<CS%MLE_density_diff) .and. (deltaRhoAtK(i)>=CS%MLE_density_diff)) then
            aFac = ( CS%MLE_density_diff - deltaRhoAtKm1(i) ) / ddRho
            MLD_fast(i,j) = dK(i) * aFac + dKm1(i) * (1. - aFac)
          endif
        enddo ! i-loop
      enddo ! k-loop
      do i = is-1, ie+1
        MLD_fast(i,j) = CS%MLE_MLD_stretch * MLD_fast(i,j)
        if ((MLD_fast(i,j)==0.) .and. (deltaRhoAtK(i)<CS%MLE_density_diff)) &
          MLD_fast(i,j) = dK(i) ! Assume mixing to the bottom
      enddo
    enddo ! j-loop
  elseif (CS%MLE_use_PBL_MLD) then
    if (.not. associated(MLD_in)) call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "Argument MLD_in was not associated!")
    do j = js-1, je+1 ; do i = is-1, ie+1
      MLD_fast(i,j) = (CS%MLE_MLD_stretch * GV%m_to_H) * MLD_in(i,j)
    enddo ; enddo
  else
    call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "No MLD to use for MLE parameterization.")
  endif

  ! Apply time filter (to remove diurnal cycle)
  if (CS%MLE_MLD_decay_time>0.) then
    if (CS%debug) then
      call hchksum(CS%MLD_filtered,'mixed_layer_restrat: MLD_filtered',G%HI,haloshift=1,scale=GV%H_to_m)
      call hchksum(MLD_in,'mixed_layer_restrat: MLD in',G%HI,haloshift=1)
    endif
    aFac = CS%MLE_MLD_decay_time / ( dt + CS%MLE_MLD_decay_time )
    bFac = dt / ( dt + CS%MLE_MLD_decay_time )
    do j = js-1, je+1 ; do i = is-1, ie+1
      ! Expression bFac*MLD_fast(i,j) + aFac*CS%MLD_filtered(i,j) is the time-filtered
      ! (running mean) of MLD. The max() allows the "running mean" to be reset
      ! instantly to a deeper MLD.
      CS%MLD_filtered(i,j) = max( MLD_fast(i,j), bFac*MLD_fast(i,j) + aFac*CS%MLD_filtered(i,j) )
      MLD_fast(i,j) = CS%MLD_filtered(i,j)
    enddo ; enddo
  endif

  ! Apply slower time filter (to remove seasonal cycle) on already filtered MLD_fast
  if (CS%MLE_MLD_decay_time2>0.) then
    if (CS%debug) then
      call hchksum(CS%MLD_filtered_slow,'mixed_layer_restrat: MLD_filtered_slow',G%HI,haloshift=1,scale=GV%H_to_m)
      call hchksum(MLD_fast,'mixed_layer_restrat: MLD fast',G%HI,haloshift=1,scale=GV%H_to_m)
    endif
    aFac = CS%MLE_MLD_decay_time2 / ( dt + CS%MLE_MLD_decay_time2 )
    bFac = dt / ( dt + CS%MLE_MLD_decay_time2 )
    do j = js-1, je+1 ; do i = is-1, ie+1
      ! Expression bFac*MLD_fast(i,j) + aFac*CS%MLD_filtered(i,j) is the time-filtered
      ! (running mean) of MLD. The max() allows the "running mean" to be reset
      ! instantly to a deeper MLD.
      CS%MLD_filtered_slow(i,j) = max( MLD_fast(i,j), bFac*MLD_fast(i,j) + aFac*CS%MLD_filtered_slow(i,j) )
      MLD_slow(i,j) = CS%MLD_filtered_slow(i,j)
    enddo ; enddo
  else
    do j = js-1, je+1 ; do i = is-1, ie+1
      MLD_slow(i,j) = MLD_fast(i,j)
    enddo ; enddo
  endif

  uDml(:) = 0.0 ; vDml(:) = 0.0
  uDml_slow(:) = 0.0 ; vDml_slow(:) = 0.0
  I4dt = 0.25 / dt
  g_Rho0 = GV%g_Earth / GV%Rho0
  h_neglect = GV%H_subroundoff
  dz_neglect = GV%H_subroundoff*GV%H_to_Z
  proper_averaging = .not. CS%MLE_use_MLD_ave_bug
  if (CS%front_length>0.) then
    res_upscale = .true.
    I_LFront = 1. / CS%front_length
  else
    res_upscale = .false.
  endif

  p0(:) = 0.0
!$OMP parallel default(none) shared(is,ie,js,je,G,GV,US,htot_fast,Rml_av_fast,tv,p0,h,h_avail,&
!$OMP                               h_neglect,g_Rho0,I4dt,CS,uhml,uhtr,dt,vhml,vhtr,   &
!$OMP                               utimescale_diag,vtimescale_diag,forces,dz_neglect, &
!$OMP                               htot_slow,MLD_slow,Rml_av_slow,VarMix,I_LFront,    &
!$OMP                               res_upscale,                                       &
!$OMP                               nz,MLD_fast,uDml_diag,vDml_diag,proper_averaging)  &
!$OMP                       private(rho_ml,h_vel,u_star,absf,mom_mixrate,timescale,    &
!$OMP                               line_is_empty, keep_going,res_scaling_fac,         &
!$OMP                               a,IhTot,b,Ihtot_slow,zpb,hAtVel,zpa,dh)            &
!$OMP                       firstprivate(uDml,vDml,uDml_slow,vDml_slow)
!$OMP do
  do j=js-1,je+1
    do i=is-1,ie+1
      htot_fast(i,j) = 0.0 ; Rml_av_fast(i,j) = 0.0
      htot_slow(i,j) = 0.0 ; Rml_av_slow(i,j) = 0.0
    enddo
    keep_going = .true.
    do k=1,nz
      do i=is-1,ie+1
        h_avail(i,j,k) = max(I4dt*G%areaT(i,j)*(h(i,j,k)-GV%Angstrom_H),0.0)
      enddo
      if (keep_going) then
        call calculate_density(tv%T(:,j,k),tv%S(:,j,k),p0,rho_ml(:),is-1,ie-is+3,tv%eqn_of_state, scale=US%kg_m3_to_R)
        line_is_empty = .true.
        do i=is-1,ie+1
          if (htot_fast(i,j) < MLD_fast(i,j)) then
            dh = h(i,j,k)
            if (proper_averaging) dh = min( h(i,j,k), MLD_fast(i,j)-htot_fast(i,j) )
            Rml_av_fast(i,j) = Rml_av_fast(i,j) + dh*rho_ml(i)
            htot_fast(i,j) = htot_fast(i,j) + dh
            line_is_empty = .false.
          endif
          if (htot_slow(i,j) < MLD_slow(i,j)) then
            dh = min( h(i,j,k), MLD_slow(i,j)-htot_slow(i,j) )
            Rml_av_slow(i,j) = Rml_av_slow(i,j) + dh*rho_ml(i)
            htot_slow(i,j) = htot_slow(i,j) + dh
            line_is_empty = .false.
          endif
        enddo
        if (line_is_empty) keep_going=.false.
      endif
    enddo

    do i=is-1,ie+1
      Rml_av_fast(i,j) = -(g_Rho0*Rml_av_fast(i,j)) / (htot_fast(i,j) + h_neglect)
      Rml_av_slow(i,j) = -(g_Rho0*Rml_av_slow(i,j)) / (htot_slow(i,j) + h_neglect)
    enddo
  enddo

  if (CS%debug) then
    call hchksum(h,'mixed_layer_restrat: h', G%HI, haloshift=1, scale=GV%H_to_m)
    call hchksum(forces%ustar,'mixed_layer_restrat: u*', G%HI, haloshift=1, scale=US%Z_to_m*US%s_to_T)
    call hchksum(MLD_fast,'mixed_layer_restrat: MLD', G%HI, haloshift=1, scale=GV%H_to_m)
    call hchksum(Rml_av_fast,'mixed_layer_restrat: rml', G%HI, haloshift=1, &
                 scale=US%m_to_Z*US%L_to_m**2*US%s_to_T**2)
  endif

! TO DO:
!   1. Mixing extends below the mixing layer to the mixed layer.  Find it!
!   2. Add exponential tail to stream-function?

!   U - Component
!$OMP do
  do j=js,je ; do I=is-1,ie
    u_star = 0.5*(forces%ustar(i,j) + forces%ustar(i+1,j))
    absf = 0.5*(abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I,J)))
    ! If needed, res_scaling_fac = min( ds, L_d ) / l_f
    if (res_upscale) res_scaling_fac = &
          ( sqrt( 0.5 * ( G%dxCu(I,j)**2 + G%dyCu(I,j)**2 ) ) * I_LFront ) &
          * min( 1., 0.5*( VarMix%Rd_dx_h(i,j) + VarMix%Rd_dx_h(i+1,j) ) )

    ! peak ML visc: u_star * 0.41 * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
    ! momentum mixing rate: pi^2*visc/h_ml^2
    ! 0.41 is the von Karmen constant, 9.8696 = pi^2.
    h_vel = 0.5*((htot_fast(i,j) + htot_fast(i+1,j)) + h_neglect) * GV%H_to_Z
    mom_mixrate = (0.41*9.8696)*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+dz_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)
    timescale = timescale * CS%ml_restrat_coef
    if (res_upscale) timescale = timescale * res_scaling_fac
    uDml(I) = timescale * G%mask2dCu(I,j)*G%dyCu(I,j)*G%IdxCu(I,j) * &
        (Rml_av_fast(i+1,j)-Rml_av_fast(i,j)) * (h_vel**2 * GV%Z_to_H)
    ! As above but using the slow filtered MLD
    h_vel = 0.5*((htot_slow(i,j) + htot_slow(i+1,j)) + h_neglect) * GV%H_to_Z
    mom_mixrate = (0.41*9.8696)*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+dz_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)
    timescale = timescale * CS%ml_restrat_coef2
    if (res_upscale) timescale = timescale * res_scaling_fac
    uDml_slow(I) = timescale * G%mask2dCu(I,j)*G%dyCu(I,j)*G%IdxCu(I,j) * &
        (Rml_av_slow(i+1,j)-Rml_av_slow(i,j)) * (h_vel**2 * GV%Z_to_H)

    if (uDml(I) + uDml_slow(I) == 0.) then
      do k=1,nz ; uhml(I,j,k) = 0.0 ; enddo
    else
      IhTot = 2.0 / ((htot_fast(i,j) + htot_fast(i+1,j)) + h_neglect)
      IhTot_slow = 2.0 / ((htot_slow(i,j) + htot_slow(i+1,j)) + h_neglect)
      zpa = 0.0 ; zpb = 0.0
      ! a(k) relates the sublayer transport to uDml with a linear profile.
      ! The sum of a(k) through the mixed layers must be 0.
      do k=1,nz
        hAtVel = 0.5*(h(i,j,k) + h(i+1,j,k))
        a(k) = PSI(zpa)                     ! Psi(z/MLD) for upper interface
        zpa = zpa - (hAtVel * IhTot)        ! z/H for lower interface
        a(k) = a(k) - PSI(zpa)              ! Transport profile
        ! Limit magnitude (uDml) if it would violate CFL
        if (a(k)*uDml(I) > 0.0) then
          if (a(k)*uDml(I) > h_avail(i,j,k)) uDml(I) = h_avail(i,j,k) / a(k)
        elseif (a(k)*uDml(I) < 0.0) then
          if (-a(k)*uDml(I) > h_avail(i+1,j,k)) uDml(I) = -h_avail(i+1,j,k) / a(k)
        endif
      enddo
      do k=1,nz
        ! Transport for slow-filtered MLD
        hAtVel = 0.5*(h(i,j,k) + h(i+1,j,k))
        b(k) = PSI(zpb)                     ! Psi(z/MLD) for upper interface
        zpb = zpb - (hAtVel * IhTot_slow)   ! z/H for lower interface
        b(k) = b(k) - PSI(zpb)              ! Transport profile
        ! Limit magnitude (uDml_slow) if it would violate CFL when added to uDml
        if (b(k)*uDml_slow(I) > 0.0) then
          if (b(k)*uDml_slow(I) > h_avail(i,j,k) - a(k)*uDml(I)) &
             uDml_slow(I) = max( 0., h_avail(i,j,k) - a(k)*uDml(I) ) / b(k)
        elseif (b(k)*uDml_slow(I) < 0.0) then
          if (-b(k)*uDml_slow(I) > h_avail(i+1,j,k) + a(k)*uDml(I)) &
             uDml_slow(I) = -max( 0., h_avail(i+1,j,k) + a(k)*uDml(I) ) / b(k)
        endif
      enddo
      do k=1,nz
        uhml(I,j,k) = a(k)*uDml(I) + b(k)*uDml_slow(I)
        uhtr(I,j,k) = uhtr(I,j,k) + uhml(I,j,k)*dt
      enddo
    endif

    utimescale_diag(I,j) = timescale
    uDml_diag(I,j) = uDml(I)
  enddo ; enddo

!  V- component
!$OMP do
  do J=js-1,je ; do i=is,ie
    u_star = 0.5*(forces%ustar(i,j) + forces%ustar(i,j+1))
    absf = 0.5*(abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))
    ! If needed, res_scaling_fac = min( ds, L_d ) / l_f
    if (res_upscale) res_scaling_fac = &
          ( sqrt( 0.5 * ( (G%dxCv(i,J))**2 + (G%dyCv(i,J))**2 ) ) * I_LFront ) &
          * min( 1., 0.5*( VarMix%Rd_dx_h(i,j) + VarMix%Rd_dx_h(i,j+1) ) )

    ! peak ML visc: u_star * 0.41 * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
    ! momentum mixing rate: pi^2*visc/h_ml^2
    ! 0.41 is the von Karmen constant, 9.8696 = pi^2.
    h_vel = 0.5*((htot_fast(i,j) + htot_fast(i,j+1)) + h_neglect) * GV%H_to_Z
    mom_mixrate = (0.41*9.8696)*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+dz_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)
    timescale = timescale * CS%ml_restrat_coef
    if (res_upscale) timescale = timescale * res_scaling_fac
    vDml(i) = timescale * G%mask2dCv(i,J)*G%dxCv(i,J)*G%IdyCv(i,J) * &
        (Rml_av_fast(i,j+1)-Rml_av_fast(i,j)) * (h_vel**2 * GV%Z_to_H)
    ! As above but using the slow filtered MLD
    h_vel = 0.5*((htot_slow(i,j) + htot_slow(i,j+1)) + h_neglect) * GV%H_to_Z
    mom_mixrate = (0.41*9.8696)*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+dz_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)
    timescale = timescale * CS%ml_restrat_coef2
    if (res_upscale) timescale = timescale * res_scaling_fac
    vDml_slow(i) = timescale * G%mask2dCv(i,J)*G%dxCv(i,J)*G%IdyCv(i,J) * &
        (Rml_av_slow(i,j+1)-Rml_av_slow(i,j)) * (h_vel**2 * GV%Z_to_H)

    if (vDml(i) + vDml_slow(i) == 0.) then
      do k=1,nz ; vhml(i,J,k) = 0.0 ; enddo
    else
      IhTot = 2.0 / ((htot_fast(i,j) + htot_fast(i,j+1)) + h_neglect)
      IhTot_slow = 2.0 / ((htot_slow(i,j) + htot_slow(i,j+1)) + h_neglect)
      zpa = 0.0 ; zpb = 0.0
      ! a(k) relates the sublayer transport to vDml with a linear profile.
      ! The sum of a(k) through the mixed layers must be 0.
      do k=1,nz
        hAtVel = 0.5*(h(i,j,k) + h(i,j+1,k))
        a(k) = PSI( zpa )                   ! Psi(z/MLD) for upper interface
        zpa = zpa - (hAtVel * IhTot)        ! z/H for lower interface
        a(k) = a(k) - PSI( zpa )            ! Transport profile
        ! Limit magnitude (vDml) if it would violate CFL
        if (a(k)*vDml(i) > 0.0) then
          if (a(k)*vDml(i) > h_avail(i,j,k)) vDml(i) = h_avail(i,j,k) / a(k)
        elseif (a(k)*vDml(i) < 0.0) then
          if (-a(k)*vDml(i) > h_avail(i,j+1,k)) vDml(i) = -h_avail(i,j+1,k) / a(k)
        endif
      enddo
      do k=1,nz
        ! Transport for slow-filtered MLD
        hAtVel = 0.5*(h(i,j,k) + h(i,j+1,k))
        b(k) = PSI(zpb)                     ! Psi(z/MLD) for upper interface
        zpb = zpb - (hAtVel * IhTot_slow)   ! z/H for lower interface
        b(k) = b(k) - PSI(zpb)              ! Transport profile
        ! Limit magnitude (vDml_slow) if it would violate CFL when added to vDml
        if (b(k)*vDml_slow(i) > 0.0) then
          if (b(k)*vDml_slow(i) > h_avail(i,j,k) - a(k)*vDml(i)) &
             vDml_slow(i) = max( 0., h_avail(i,j,k) - a(k)*vDml(i) ) / b(k)
        elseif (b(k)*vDml_slow(i) < 0.0) then
          if (-b(k)*vDml_slow(i) > h_avail(i,j+1,k) + a(k)*vDml(i)) &
             vDml_slow(i) = -max( 0., h_avail(i,j+1,k) + a(k)*vDml(i) ) / b(k)
        endif
      enddo
      do k=1,nz
        vhml(i,J,k) = a(k)*vDml(i) + b(k)*vDml_slow(i)
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
  if (CS%id_uhml > 0 .or. CS%id_vhml > 0) &
    ! Remapped uhml and vhml require east/north halo updates of h
    call pass_var(h, G%domain, To_West+To_South+Omit_Corners, halo=1)
  call diag_update_remap_grids(CS%diag)

  ! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_urestrat_time > 0) call post_data(CS%id_urestrat_time, utimescale_diag, CS%diag)
    if (CS%id_vrestrat_time > 0) call post_data(CS%id_vrestrat_time, vtimescale_diag, CS%diag)
    if (CS%id_uhml          > 0) call post_data(CS%id_uhml, uhml, CS%diag)
    if (CS%id_vhml          > 0) call post_data(CS%id_vhml, vhml, CS%diag)
    if (CS%id_MLD           > 0) call post_data(CS%id_MLD, MLD_fast, CS%diag)
    if (CS%id_Rml           > 0) call post_data(CS%id_Rml, Rml_av_fast, CS%diag)
    if (CS%id_uDml          > 0) call post_data(CS%id_uDml, uDml_diag, CS%diag)
    if (CS%id_vDml          > 0) call post_data(CS%id_vDml, vDml_diag, CS%diag)

    if (CS%id_uml > 0) then
      do J=js,je ; do i=is-1,ie
        h_vel = 0.5*((htot_fast(i,j) + htot_fast(i+1,j)) + h_neglect)
        uDml_diag(I,j) = uDml_diag(I,j) / (0.01*h_vel) * G%IdyCu(I,j) * (PSI(0.)-PSI(-.01))
      enddo ; enddo
      call post_data(CS%id_uml, uDml_diag, CS%diag)
    endif
    if (CS%id_vml > 0) then
      do J=js-1,je ; do i=is,ie
        h_vel = 0.5*((htot_fast(i,j) + htot_fast(i,j+1)) + h_neglect)
        vDml_diag(i,J) = vDml_diag(i,J) / (0.01*h_vel) * G%IdxCv(i,J) * (PSI(0.)-PSI(-.01))
      enddo ; enddo
      call post_data(CS%id_vml, vDml_diag, CS%diag)
    endif
  endif
  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  ! This needs to happen after the H update and before the next post_data.
  call diag_update_remap_grids(CS%diag)

end subroutine mixedlayer_restrat_general


!> Calculates a restratifying flow assuming a 2-layer bulk mixed layer.
subroutine mixedlayer_restrat_BML(h, uhtr, vhtr, tv, forces, dt, G, GV, US, CS)
  type(ocean_grid_type),                     intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                     !!   [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                     !!   [H L2 ~> m3 or kg]
  type(thermo_var_ptrs),                     intent(in)    :: tv     !< Thermodynamic variables structure
  type(mech_forcing),                        intent(in)    :: forces !< A structure with the driving mechanical forces
  real,                                      intent(in)    :: dt     !< Time increment [T ~> s]
  type(mixedlayer_restrat_CS),               pointer       :: CS     !< Module control structure
  ! Local variables
  real :: uhml(SZIB_(G),SZJ_(G),SZK_(G)) ! zonal mixed layer transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: vhml(SZI_(G),SZJB_(G),SZK_(G)) ! merid mixed layer transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    h_avail               ! The volume available for diffusion out of each face of each
                          ! sublayer of the mixed layer, divided by dt [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    htot, &               ! The sum of the thicknesses of layers in the mixed layer [H ~> m or kg m-2]
    Rml_av                ! g_Rho0 times the average mixed layer density [L2 Z-1 T-2 ~> m s-2]
  real :: g_Rho0          ! G_Earth/Rho0 [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1]
  real :: Rho0(SZI_(G))   ! Potential density relative to the surface [R ~> kg m-3]
  real :: p0(SZI_(G))     ! A pressure of 0 [Pa]

  real :: h_vel           ! htot interpolated onto velocity points [Z ~> m]. (The units are not H.)
  real :: absf            ! absolute value of f, interpolated to velocity points [T-1 ~> s-1]
  real :: u_star          ! surface friction velocity, interpolated to velocity points [Z T-1 ~> m s-1].
  real :: mom_mixrate     ! rate at which momentum is homogenized within mixed layer [T-1 ~> s-1]
  real :: timescale       ! mixing growth timescale [T ~> s]
  real :: h_neglect       ! tiny thickness usually lost in roundoff and can be neglected [H ~> m or kg m-2]
  real :: dz_neglect      ! tiny thickness that usually lost in roundoff and can be neglected [Z ~> m]
  real :: I4dt            ! 1/(4 dt) [T-1 ~> s-1]
  real :: I2htot          ! Twice the total mixed layer thickness at velocity points [H ~> m or kg m-2]
  real :: z_topx2         ! depth of the top of a layer at velocity points [H ~> m or kg m-2]
  real :: hx2             ! layer thickness at velocity points [H ~> m or kg m-2]
  real :: a(SZK_(G))      ! A non-dimensional value relating the overall flux
                          ! magnitudes (uDml & vDml) to the realized flux in a
                          ! layer.  The vertical sum of a() through the pieces of
                          ! the mixed layer must be 0.
  real :: uDml(SZIB_(G))  ! The zonal and meridional volume fluxes in the upper
  real :: vDml(SZI_(G))   ! half of the mixed layer [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: utimescale_diag(SZIB_(G),SZJ_(G)) ! The restratification timescales
  real :: vtimescale_diag(SZI_(G),SZJB_(G)) ! in the zonal and meridional
                                            ! directions [T ~> s], stored in 2-D
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
  g_Rho0     = GV%g_Earth / GV%Rho0
  use_EOS    = associated(tv%eqn_of_state)
  h_neglect  = GV%H_subroundoff
  dz_neglect = GV%H_subroundoff*GV%H_to_Z

  if (.not.use_EOS) call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "An equation of state must be used with this module.")

  ! Fix this later for nkml >= 3.

  p0(:) = 0.0
!$OMP parallel default(none) shared(is,ie,js,je,G,GV,US,htot,Rml_av,tv,p0,h,h_avail,   &
!$OMP                               h_neglect,g_Rho0,I4dt,CS,uhml,uhtr,dt,vhml,vhtr,   &
!$OMP                               utimescale_diag,vtimescale_diag,forces,dz_neglect, &
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
      call calculate_density(tv%T(:,j,k),tv%S(:,j,k),p0,Rho0(:),is-1,ie-is+3,tv%eqn_of_state, scale=US%kg_m3_to_R)
      do i=is-1,ie+1
        Rml_av(i,j) = Rml_av(i,j) + h(i,j,k)*Rho0(i)
        htot(i,j) = htot(i,j) + h(i,j,k)
        h_avail(i,j,k) = max(I4dt*G%areaT(i,j)*(h(i,j,k)-GV%Angstrom_H),0.0)
      enddo
    enddo

    do i=is-1,ie+1
      Rml_av(i,j) = (g_Rho0*Rml_av(i,j)) / (htot(i,j) + h_neglect)
    enddo
  enddo

! TO DO:
!   1. Mixing extends below the mixing layer to the mixed layer.  Find it!
!   2. Add exponential tail to stream-function?

!   U - Component
!$OMP do
  do j=js,je; do I=is-1,ie
    h_vel = 0.5*(htot(i,j) + htot(i+1,j)) * GV%H_to_Z

    u_star = 0.5*(forces%ustar(i,j) + forces%ustar(i+1,j))
    absf = 0.5*(abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I,J)))
    ! peak ML visc: u_star * 0.41 * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
    ! momentum mixing rate: pi^2*visc/h_ml^2
    ! 0.41 is the von Karmen constant, 9.8696 = pi^2.
    mom_mixrate = (0.41*9.8696)*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+dz_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)

    timescale = timescale * CS%ml_restrat_coef
!      timescale = timescale*(2?)*(L_def/L_MLI) * min(EKE/MKE,1.0 + (G%dyCv(i,j)/L_def)**2)

    uDml(I) = timescale * G%mask2dCu(I,j)*G%dyCu(I,j)*G%IdxCu(I,j) * &
        (Rml_av(i+1,j)-Rml_av(i,j)) * (h_vel**2 * GV%Z_to_H)

    if (uDml(I) == 0) then
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

    uDml_diag(I,j) = uDml(I)
    utimescale_diag(I,j) = timescale
  enddo; enddo

!  V- component
!$OMP do
  do J=js-1,je ; do i=is,ie
    h_vel = 0.5*(htot(i,j) + htot(i,j+1)) * GV%H_to_Z

    u_star = 0.5*(forces%ustar(i,j) + forces%ustar(i,j+1))
    absf = 0.5*(abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))
    ! peak ML visc: u_star * 0.41 * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
    ! momentum mixing rate: pi^2*visc/h_ml^2
    ! 0.41 is the von Karmen constant, 9.8696 = pi^2.
    mom_mixrate = (0.41*9.8696)*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+dz_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)

    timescale = timescale * CS%ml_restrat_coef
!     timescale = timescale*(2?)*(L_def/L_MLI) * min(EKE/MKE,1.0 + (G%dyCv(i,j)/L_def)**2)

    vDml(i) = timescale * G%mask2dCv(i,J)*G%dxCv(i,J)*G%IdyCv(i,J) * &
        (Rml_av(i,j+1)-Rml_av(i,j)) * (h_vel**2 * GV%Z_to_H)
    if (vDml(i) == 0) then
      do k=1,nkml ; vhml(i,J,k) = 0.0 ; enddo
    else
      I2htot = 1.0 / (htot(i,j) + htot(i,j+1) + h_neglect)
      z_topx2 = 0.0
      ! a(k) relates the sublayer transport to vDml with a linear profile.
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

    vtimescale_diag(i,J) = timescale
    vDml_diag(i,J) = vDml(i)
  enddo; enddo

!$OMP do
  do j=js,je ; do k=1,nkml ; do i=is,ie
    h(i,j,k) = h(i,j,k) - dt*G%IareaT(i,j) * &
        ((uhml(I,j,k) - uhml(I-1,j,k)) + (vhml(i,J,k) - vhml(i,J-1,k)))
  enddo ; enddo ; enddo
!$OMP end parallel

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  if (CS%id_uhml > 0 .or. CS%id_vhml > 0) &
    ! Remapped uhml and vhml require east/north halo updates of h
    call pass_var(h, G%domain, To_West+To_South+Omit_Corners, halo=1)
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


!> Initialize the mixed layer restratification module
logical function mixedlayer_restrat_init(Time, G, GV, US, param_file, diag, CS, restart_CS)
  type(time_type),             intent(in)    :: Time       !< Current model time
  type(ocean_grid_type),       intent(inout) :: G          !< Ocean grid structure
  type(verticalGrid_type),     intent(in)    :: GV         !< Ocean vertical grid structure
  type(unit_scale_type),       intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),       intent(in)    :: param_file !< Parameter file to parse
  type(diag_ctrl), target,     intent(inout) :: diag       !< Regulate diagnostics
  type(mixedlayer_restrat_CS), pointer       :: CS         !< Module control structure
  type(MOM_restart_CS),        pointer       :: restart_CS !< A pointer to the restart control structure

  ! Local variables
  real :: H_rescale  ! A rescaling factor for thicknesses from the representation in
                     ! a restart file to the internal representation in this run.
  real :: flux_to_kg_per_s ! A unit conversion factor for fluxes.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  integer :: i, j

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MIXEDLAYER_RESTRAT", mixedlayer_restrat_init, &
             "If true, a density-gradient dependent re-stratifying "//&
             "flow is imposed in the mixed layer. Can be used in ALE mode "//&
             "without restriction but in layer mode can only be used if "//&
             "BULKMIXEDLAYER is true.", default=.false.)
  if (.not. mixedlayer_restrat_init) return

  if (.not.associated(CS)) then
    call MOM_error(FATAL, "mixedlayer_restrat_init called without an associated control structure.")
  endif

  ! Nonsense values to cause problems when these parameters are not used
  CS%MLE_MLD_decay_time = -9.e9*US%s_to_T
  CS%MLE_density_diff = -9.e9*US%kg_m3_to_R
  CS%MLE_tail_dh = -9.e9
  CS%MLE_use_PBL_MLD = .false.
  CS%MLE_MLD_stretch = -9.e9

  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "FOX_KEMPER_ML_RESTRAT_COEF", CS%ml_restrat_coef, &
             "A nondimensional coefficient that is proportional to "//&
             "the ratio of the deformation radius to the dominant "//&
             "lengthscale of the submesoscale mixed layer "//&
             "instabilities, times the minimum of the ratio of the "//&
             "mesoscale eddy kinetic energy to the large-scale "//&
             "geostrophic kinetic energy or 1 plus the square of the "//&
             "grid spacing over the deformation radius, as detailed "//&
             "by Fox-Kemper et al. (2010)", units="nondim", default=0.0)
  ! We use GV%nkml to distinguish between the old and new implementation of MLE.
  ! The old implementation only works for the layer model with nkml>0.
  if (GV%nkml==0) then
    call get_param(param_file, mdl, "FOX_KEMPER_ML_RESTRAT_COEF2", CS%ml_restrat_coef2, &
             "As for FOX_KEMPER_ML_RESTRAT_COEF but used in a second application "//&
             "of the MLE restratification parameterization.", units="nondim", default=0.0)
    call get_param(param_file, mdl, "MLE_FRONT_LENGTH", CS%front_length, &
             "If non-zero, is the frontal-length scale used to calculate the "//&
             "upscaling of buoyancy gradients that is otherwise represented "//&
             "by the parameter FOX_KEMPER_ML_RESTRAT_COEF. If MLE_FRONT_LENGTH is "//&
             "non-zero, it is recommended to set FOX_KEMPER_ML_RESTRAT_COEF=1.0.",&
             units="m", default=0.0, scale=US%m_to_L)
    call get_param(param_file, mdl, "MLE_USE_PBL_MLD", CS%MLE_use_PBL_MLD, &
             "If true, the MLE parameterization will use the mixed-layer "//&
             "depth provided by the active PBL parameterization. If false, "//&
             "MLE will estimate a MLD based on a density difference with the "//&
             "surface using the parameter MLE_DENSITY_DIFF.", default=.false.)
    call get_param(param_file, mdl, "MLE_MLD_DECAY_TIME", CS%MLE_MLD_decay_time, &
             "The time-scale for a running-mean filter applied to the mixed-layer "//&
             "depth used in the MLE restratification parameterization. When "//&
             "the MLD deepens below the current running-mean the running-mean "//&
             "is instantaneously set to the current MLD.", units="s", default=0., scale=US%s_to_T)
    call get_param(param_file, mdl, "MLE_MLD_DECAY_TIME2", CS%MLE_MLD_decay_time2, &
             "The time-scale for a running-mean filter applied to the filtered "//&
             "mixed-layer depth used in a second MLE restratification parameterization. "//&
             "When the MLD deepens below the current running-mean the running-mean "//&
             "is instantaneously set to the current MLD.", units="s", default=0., scale=US%s_to_T)
    if (.not. CS%MLE_use_PBL_MLD) then
      call get_param(param_file, mdl, "MLE_DENSITY_DIFF", CS%MLE_density_diff, &
             "Density difference used to detect the mixed-layer "//&
             "depth used for the mixed-layer eddy parameterization "//&
             "by Fox-Kemper et al. (2010)", units="kg/m3", default=0.03, scale=US%kg_m3_to_R)
    endif
    call get_param(param_file, mdl, "MLE_TAIL_DH", CS%MLE_tail_dh, &
             "Fraction by which to extend the mixed-layer restratification "//&
             "depth used for a smoother stream function at the base of "//&
             "the mixed-layer.", units="nondim", default=0.0)
    call get_param(param_file, mdl, "MLE_MLD_STRETCH", CS%MLE_MLD_stretch, &
             "A scaling coefficient for stretching/shrinking the MLD "//&
             "used in the MLE scheme. This simply multiplies MLD wherever used.",&
             units="nondim", default=1.0)
    call get_param(param_file, mdl, "MLE_USE_MLD_AVE_BUG", CS%MLE_use_MLD_ave_bug, &
             "If true, do not account for MLD mismatch to interface positions.",&
             default=.false.)
  endif

  CS%diag => diag

  flux_to_kg_per_s = GV%H_to_kg_m2 * US%L_to_m**2 * US%s_to_T

  CS%id_uhml = register_diag_field('ocean_model', 'uhml', diag%axesCuL, Time, &
      'Zonal Thickness Flux to Restratify Mixed Layer', 'kg s-1', &
      conversion=flux_to_kg_per_s, y_cell_method='sum', v_extensive=.true.)
  CS%id_vhml = register_diag_field('ocean_model', 'vhml', diag%axesCvL, Time, &
      'Meridional Thickness Flux to Restratify Mixed Layer', 'kg s-1', &
      conversion=flux_to_kg_per_s, x_cell_method='sum', v_extensive=.true.)
  CS%id_urestrat_time = register_diag_field('ocean_model', 'MLu_restrat_time', diag%axesCu1, Time, &
      'Mixed Layer Zonal Restratification Timescale', 's', conversion=US%T_to_s)
  CS%id_vrestrat_time = register_diag_field('ocean_model', 'MLv_restrat_time', diag%axesCv1, Time, &
      'Mixed Layer Meridional Restratification Timescale', 's', conversion=US%T_to_s)
  CS%id_MLD = register_diag_field('ocean_model', 'MLD_restrat', diag%axesT1, Time, &
      'Mixed Layer Depth as used in the mixed-layer restratification parameterization', 'm', &
      conversion=GV%H_to_m)
  CS%id_Rml = register_diag_field('ocean_model', 'ML_buoy_restrat', diag%axesT1, Time, &
      'Mixed Layer Buoyancy as used in the mixed-layer restratification parameterization', &
      'm s2', conversion=US%m_to_Z*(US%L_to_m**2)*(US%s_to_T**2))
  CS%id_uDml = register_diag_field('ocean_model', 'udml_restrat', diag%axesCu1, Time, &
      'Transport stream function amplitude for zonal restratification of mixed layer', &
      'm3 s-1', conversion=GV%H_to_m*(US%L_to_m**2)*US%s_to_T)
  CS%id_vDml = register_diag_field('ocean_model', 'vdml_restrat', diag%axesCv1, Time, &
      'Transport stream function amplitude for meridional restratification of mixed layer', &
      'm3 s-1', conversion=GV%H_to_m*(US%L_to_m**2)*US%s_to_T)
  CS%id_uml = register_diag_field('ocean_model', 'uml_restrat', diag%axesCu1, Time, &
      'Surface zonal velocity component of mixed layer restratification', &
      'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vml = register_diag_field('ocean_model', 'vml_restrat', diag%axesCv1, Time, &
      'Surface meridional velocity component of mixed layer restratification', &
      'm s-1', conversion=US%L_T_to_m_s)

  ! Rescale variables from restart files if the internal dimensional scalings have changed.
  if (CS%MLE_MLD_decay_time>0. .or. CS%MLE_MLD_decay_time2>0.) then
    if (query_initialized(CS%MLD_filtered, "MLD_MLE_filtered", restart_CS) .and. &
        (GV%m_to_H_restart /= 0.0) .and. (GV%m_to_H_restart /= GV%m_to_H)) then
      H_rescale = GV%m_to_H / GV%m_to_H_restart
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        CS%MLD_filtered(i,j) = H_rescale * CS%MLD_filtered(i,j)
      enddo ; enddo
    endif
  endif
  if (CS%MLE_MLD_decay_time2>0.) then
    if (query_initialized(CS%MLD_filtered_slow, "MLD_MLE_filtered_slow", restart_CS) .and. &
        (GV%m_to_H_restart /= 0.0) .and. (GV%m_to_H_restart /= GV%m_to_H)) then
      H_rescale = GV%m_to_H / GV%m_to_H_restart
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        CS%MLD_filtered_slow(i,j) = H_rescale * CS%MLD_filtered_slow(i,j)
      enddo ; enddo
    endif
  endif

  ! If MLD_filtered is being used, we need to update halo regions after a restart
  if (associated(CS%MLD_filtered)) call pass_var(CS%MLD_filtered, G%domain)

end function mixedlayer_restrat_init

!> Allocate and register fields in the mixed layer restratification structure for restarts
subroutine mixedlayer_restrat_register_restarts(HI, param_file, CS, restart_CS)
  ! Arguments
  type(hor_index_type),        intent(in)    :: HI         !< Horizontal index structure
  type(param_file_type),       intent(in)    :: param_file !< Parameter file to parse
  type(mixedlayer_restrat_CS), pointer       :: CS         !< Module control structure
  type(MOM_restart_CS),        pointer       :: restart_CS !< A pointer to the restart control structure
  ! Local variables
  type(vardesc) :: vd
  logical :: mixedlayer_restrat_init

  ! Check to see if this module will be used
  call get_param(param_file, mdl, "MIXEDLAYER_RESTRAT", mixedlayer_restrat_init, &
             default=.false., do_not_log=.true.)
  if (.not. mixedlayer_restrat_init) return

  ! Allocate the control structure. CS will be later populated by mixedlayer_restrat_init()
  if (associated(CS)) call MOM_error(FATAL, &
       "mixedlayer_restrat_register_restarts called with an associated control structure.")
  allocate(CS)

  call get_param(param_file, mdl, "MLE_MLD_DECAY_TIME", CS%MLE_MLD_decay_time, &
                 default=0., do_not_log=.true.)
  call get_param(param_file, mdl, "MLE_MLD_DECAY_TIME2", CS%MLE_MLD_decay_time2, &
                 default=0., do_not_log=.true.)
  if (CS%MLE_MLD_decay_time>0. .or. CS%MLE_MLD_decay_time2>0.) then
    ! CS%MLD_filtered is used to keep a running mean of the PBL's actively mixed MLD.
    allocate(CS%MLD_filtered(HI%isd:HI%ied,HI%jsd:HI%jed)) ; CS%MLD_filtered(:,:) = 0.
    vd = var_desc("MLD_MLE_filtered","m","Time-filtered MLD for use in MLE", &
                  hor_grid='h', z_grid='1')
    call register_restart_field(CS%MLD_filtered, vd, .false., restart_CS)
  endif
  if (CS%MLE_MLD_decay_time2>0.) then
    ! CS%MLD_filtered_slow is used to keep a running mean of the PBL's seasonal or winter MLD.
    allocate(CS%MLD_filtered_slow(HI%isd:HI%ied,HI%jsd:HI%jed)) ; CS%MLD_filtered_slow(:,:) = 0.
    vd = var_desc("MLD_MLE_filtered_slow","m","c Slower time-filtered MLD for use in MLE", &
                  hor_grid='h', z_grid='1')
    call register_restart_field(CS%MLD_filtered_slow, vd, .false., restart_CS)
  endif

end subroutine mixedlayer_restrat_register_restarts

!> \namespace mom_mixed_layer_restrat
!!
!! \section section_mle Mixed-layer eddy parameterization module
!!
!! The subroutines in this file implement a parameterization of unresolved viscous
!! mixed layer restratification of the mixed layer as described in Fox-Kemper et
!! al., 2008, and whose impacts are described in Fox-Kemper et al., 2011.
!! This is derived in part from the older parameterization that is described in
!! Hallberg (Aha Hulikoa, 2003), which this new parameterization surpasses, which
!! in turn is based on the sub-inertial mixed layer theory of Young (JPO, 1994).
!! There is no net horizontal volume transport due to this parameterization, and
!! no direct effect below the mixed layer.
!!
!! This parameterization sets the restratification timescale to agree with
!! high-resolution studies of mixed layer restratification.
!!
!! The run-time parameter FOX_KEMPER_ML_RESTRAT_COEF is a non-dimensional number of
!! order a few tens, proportional to the ratio of the deformation radius or the
!! grid scale (whichever is smaller to the dominant horizontal length-scale of the
!! sub-meso-scale mixed layer instabilities.
!!
!! \subsection section_mle_nutshell "Sub-meso" in a nutshell
!!
!! The parameterization is colloquially referred to as "sub-meso".
!!
!! The original Fox-Kemper et al., (2008b) paper proposed a quasi-Stokes
!! advection described by the stream function (eq. 5 of Fox-Kemper et al., 2011):
!! \f[
!!    {\bf \Psi}_o = C_e \frac{ H^2 \nabla \bar{b} \times \hat{\bf z} }{ |f| } \mu(z)
!! \f]
!!
!! where the vertical profile function is
!! \f[
!!    \mu(z) = \max \left\{ 0, \left[ 1 - \left(\frac{2z}{H}+1\right)^2 \right]
!!                            \left[ 1 + \frac{5}{21} \left(\frac{2z}{H}+1\right)^2 \right] \right\}
!! \f]
!! and \f$ H \f$ is the mixed-layer depth, \f$ f \f$ is the local Coriolis parameter, \f$ C_e \sim 0.06-0.08 \f$ and
!! \f$ \nabla \bar{b} \f$ is a depth mean buoyancy gradient averaged over the mixed layer.
!!
!! For use in coarse-resolution models, an upscaling of the buoyancy gradients and adaption for the equator
!! leads to the following parameterization (eq. 6 of Fox-Kemper et al., 2011):
!! \f[
!!    {\bf \Psi} = C_e \Gamma_\Delta \frac{\Delta s}{l_f} \frac{ H^2 \nabla \bar{b} \times \hat{\bf z} }
!!                 { \sqrt{ f^2 + \tau^{-2}} } \mu(z)
!! \f]
!! where \f$ \Delta s \f$ is the minimum of grid-scale and deformation radius,
!! \f$ l_f \f$ is the width of the mixed-layer fronts, and \f$ \Gamma_\Delta=1 \f$.
!! \f$ \tau \f$ is a time-scale for mixing momentum across the mixed layer.
!! \f$ l_f \f$ is thought to be of order hundreds of meters.
!!
!! The upscaling factor \f$ \frac{\Delta s}{l_f} \f$ can be a global constant, model parameter FOX_KEMPER_ML_RESTRAT,
!! so that in practice the parameterization is:
!! \f[
!!    {\bf \Psi} = C_e \Gamma_\Delta \frac{ H^2 \nabla \bar{b} \times \hat{\bf z} }{ \sqrt{ f^2 + \tau^{-2}} } \mu(z)
!! \f]
!! with non-unity \f$ \Gamma_\Delta \f$.
!!
!! \f$ C_e \f$ is hard-coded as 0.0625. \f$ \tau \f$ is calculated from the surface friction velocity \f$ u^* \f$.
!! \todo Explain expression for momentum mixing time-scale.
!!
!! \subsection section_mle_filtering Time-filtering of mixed-layer depth
!!
!! Using the instantaneous mixed-layer depth is inconsistent with the finite life-time of
!! mixed-layer instabilities. We provide a one-sided running-mean filter of mixed-layer depth, \f$ H \f$, of the form:
!! \f[
!!    \bar{H} \leftarrow \max \left( H, \frac{ \Delta t H + \tau_h \bar{H} }{ \Delta t + \tau_h } \right)
!! \f]
!! which allows the effective mixed-layer depth seen by the parameterization, \f$\bar{H}\f$, to instantaneously deepen
!! but to decay with time-scale \f$ \tau_h \f$.
!! \f$ \bar{H} \f$ is substituted for \f$ H \f$ in the above equations.
!!
!! \subsection section_mle_mld Defining the mixed-layer-depth
!!
!! If the parameter MLE_USE_PBL_MLD=True then the mixed-layer depth is defined/diagnosed by the
!! boundary-layer parameterization (e.g. ePBL, KPP, etc.).
!!
!! If the parameter MLE_USE_PBL_MLD=False then the mixed-layer depth is diagnosed in this module
!! as the depth of a given density difference, \f$ \Delta \rho \f$, with the surface where the
!! density difference is the parameter MLE_DENSITY_DIFF.
!!
!! \subsection section_mle_ref References
!!
!! Fox-Kemper, B., Ferrari, R. and Hallberg, R., 2008:
!! Parameterization of Mixed Layer Eddies. Part I: Theory and Diagnosis
!! J. Phys. Oceangraphy, 38 (6), p1145-1165.
!! https://doi.org/10.1175/2007JPO3792.1
!!
!! Fox-Kemper, B. and Ferrari, R. 2008:
!! Parameterization of Mixed Layer Eddies. Part II: Prognosis and Impact
!! J. Phys. Oceangraphy, 38 (6), p1166-1179.
!! https://doi.org/10.1175/2007JPO3788.1
!!
!! B. Fox-Kemper, G. Danabasoglu, R. Ferrari, S.M. Griffies, R.W. Hallberg, M.M. Holland, M.E. Maltrud,
!! S. Peacock, and B.L. Samuels, 2011: Parameterization of mixed layer eddies. III: Implementation and impact
!! in global ocean climate simulations. Ocean Modell., 39(1), p61-78.
!! https://doi.org/10.1016/j.ocemod.2010.09.002
!!
!! | Symbol                       | Module parameter      |
!! | ---------------------------- | --------------------- |
!! | \f$ \Gamma_\Delta \f$        | FOX_KEMPER_ML_RESTRAT |
!! | \f$ l_f \f$                  | MLE_FRONT_LENGTH      |
!! | \f$ \tau_h \f$               | MLE_MLD_DECAY_TIME    |
!! | \f$ \Delta \rho \f$          | MLE_DENSITY_DIFF      |

end module MOM_mixed_layer_restrat
