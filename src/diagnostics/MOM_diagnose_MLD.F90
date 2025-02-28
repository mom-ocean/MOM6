!> Provides functions for some diabatic processes such as fraxil, brine rejection,
!! tendency due to surface flux divergence.
module MOM_diagnose_mld

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data
use MOM_diag_mediator, only : diag_ctrl
use MOM_EOS,           only : calculate_density, calculate_TFreeze, EOS_domain
use MOM_EOS,           only : calculate_specific_vol_derivs, calculate_density_derivs
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_grid,          only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public diagnoseMLDbyEnergy, diagnoseMLDbyDensityDifference

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains
!> Diagnose a mixed layer depth (MLD) determined by a given density difference with the surface.
!> This routine is appropriate in MOM_diabatic_aux due to its position within the time stepping.
subroutine diagnoseMLDbyDensityDifference(id_MLD, h, tv, densityDiff, G, GV, US, diagPtr, &
                                          ref_h_mld, id_ref_z, id_ref_rho, id_N2subML, id_MLDsq, dz_subML)
  type(ocean_grid_type),   intent(in) :: G           !< Grid type
  type(verticalGrid_type), intent(in) :: GV          !< ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US          !< A dimensional unit scaling type
  integer,                 intent(in) :: id_MLD      !< Handle (ID) of MLD diagnostic
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h           !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in) :: tv          !< Structure containing pointers to any
                                                     !! available thermodynamic fields.
  real,                    intent(in) :: densityDiff !< Density difference to determine MLD [R ~> kg m-3]
  type(diag_ctrl),         pointer    :: diagPtr     !< Diagnostics structure
  real,                    intent(in) :: ref_h_mld   !< Depth of the calculated "surface" densisty [Z ~> m]
  integer,                 intent(in) :: id_ref_z    !< Handle (ID) of reference depth diagnostic
  integer,                 intent(in) :: id_ref_rho  !< Handle (ID) of reference density diagnostic
  integer,       optional, intent(in) :: id_N2subML  !< Optional handle (ID) of subML stratification
  integer,       optional, intent(in) :: id_MLDsq    !< Optional handle (ID) of squared MLD
  real,          optional, intent(in) :: dz_subML    !< The distance over which to calculate N2subML
                                                     !! or 50 m if missing [Z ~> m]

  ! Local variables
  real, dimension(SZI_(G)) :: deltaRhoAtKm1, deltaRhoAtK ! Density differences [R ~> kg m-3].
  real, dimension(SZI_(G)) :: pRef_MLD, pRef_N2 ! Reference pressures [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G)) :: hRef_MLD          ! Reference depth  [Z ~> m].
  real, dimension(SZI_(G)) :: H_subML, dH_N2    ! Summed thicknesses used in N2 calculation [H ~> m or kg m-2]
  real, dimension(SZI_(G)) :: dZ_N2             ! Summed vertical distance used in N2 calculation [Z ~> m]
  real, dimension(SZI_(G)) :: T_subML, T_deeper ! Temperatures used in the N2 calculation [C ~> degC].
  real, dimension(SZI_(G)) :: S_subML, S_deeper ! Salinities used in the N2 calculation [S ~> ppt].
  real, dimension(SZI_(G)) :: rho_subML, rho_deeper ! Densities used in the N2 calculation [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)) :: dZ_2d   ! Layer thicknesses in depth units [Z ~> m]
  real, dimension(SZI_(G)) :: dZ, dZm1         ! Layer thicknesses associated with interfaces [Z ~> m]
  real, dimension(SZI_(G)) :: rhoSurf          ! Density used in finding the mixed layer depth [R ~> kg m-3].
  real, dimension(SZI_(G), SZJ_(G)) :: z_ref_diag ! The actual depth of the reference density [Z ~> m].
  real, dimension(SZI_(G), SZJ_(G)) :: MLD     ! Diagnosed mixed layer depth [Z ~> m].
  real, dimension(SZI_(G), SZJ_(G)) :: subMLN2 ! Diagnosed stratification below ML [T-2 ~> s-2].
  real, dimension(SZI_(G), SZJ_(G)) :: MLD2    ! Diagnosed MLD^2 [Z2 ~> m2].
  logical, dimension(SZI_(G)) :: N2_region_set ! If true, all necessary values for calculating N2
                                               ! have been stored already.
  real :: gE_Rho0          ! The gravitational acceleration, sometimes divided by the Boussinesq
                           ! reference density [H T-2 R-1 ~> m4 s-2 kg-1 or m s-2].
  real :: dZ_sub_ML        ! Depth below ML over which to diagnose stratification [Z ~> m]
  real :: aFac             ! A nondimensional factor [nondim]
  real :: ddRho            ! A density difference [R ~> kg m-3]
  real :: dddpth           ! A depth difference [Z ~> m]
  real :: rhoSurf_k, rhoSurf_km1  ! Desisty in the layers below and above the target reference depth [R ~> kg m-3].
  real, dimension(SZI_(G), SZJ_(G)) :: rhoSurf_2d ! The density that is considered the "surface" when calculating
                                                  ! the MLD. It can be saved as a diagnostic [R ~> kg m-3].
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, is, ie, js, je, k, nz, id_N2, id_SQ

  id_SQ = -1 ; if (PRESENT(id_MLDsq)) id_SQ = id_MLDsq

  id_N2 = -1
  if (present(id_N2subML)) then
    if (present(dz_subML)) then
      id_N2 = id_N2subML
      dZ_sub_ML = dz_subML
    else
      call MOM_error(FATAL, "When the diagnostic of the subML stratification is "//&
                "requested by providing id_N2_subML to diagnoseMLDbyDensityDifference, "//&
                "the distance over which to calculate that distance must also be provided.")
    endif
  endif

  gE_rho0 = (US%L_to_Z**2 * GV%g_Earth) / GV%H_to_RZ

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  hRef_MLD(:) = ref_h_mld
  pRef_MLD(:) = GV%H_to_RZ*GV%g_Earth*ref_h_mld
  z_ref_diag(:,:) = 0.

  EOSdom(:) = EOS_domain(G%HI)
  do j=js,je
    ! Find the vertical distances across layers.
    call thickness_to_dz(h, tv, dZ_2d, j, G, GV)

    if (pRef_MLD(is) /= 0.0) then
      rhoSurf(:) = 0.0
      do i=is,ie
        dZ(i) = 0.5 * dZ_2d(i,1) ! Depth of center of surface layer
        if (dZ(i) >= hRef_MLD(i)) then
          call calculate_density(tv%T(i,j,1), tv%S(i,j,1), pRef_MLD(i), rhoSurf_k, tv%eqn_of_state)
          rhoSurf(i) = rhoSurf_k
        endif
      enddo
      do k=2,nz
        do i=is,ie
          dZm1(i) = dZ(i) ! Depth of center of layer K-1
          dZ(i) = dZ(i) + 0.5 * ( dZ_2d(i,k) + dZ_2d(i,k-1) ) ! Depth of center of layer K
          dddpth = dZ(i) - dZm1(i)
          if ((rhoSurf(i) == 0.) .and. &
              (dZm1(i) < hRef_MLD(i)) .and. (dZ(i) >= hRef_MLD(i))) then
            aFac = ( hRef_MLD(i) - dZm1(i) ) / dddpth
            z_ref_diag(i,j) = (dZ(i) * aFac + dZm1(i) * (1. - aFac))
            call calculate_density(tv%T(i,j,k)  , tv%S(i,j,k)  , pRef_MLD(i), rhoSurf_k, tv%eqn_of_state)
            call calculate_density(tv%T(i,j,k-1), tv%S(i,j,k-1), pRef_MLD(i), rhoSurf_km1, tv%eqn_of_state)
            rhoSurf(i) = (rhoSurf_k * aFac + rhoSurf_km1 * (1. - aFac))
            H_subML(i) = h(i,j,k)
          elseif ((rhoSurf(i) == 0.) .and. (k >= nz)) then
            call calculate_density(tv%T(i,j,1), tv%S(i,j,1), pRef_MLD(i), rhoSurf_k, tv%eqn_of_state)
            rhoSurf(i) = rhoSurf_k
          endif
        enddo
      enddo
      do i=is,ie
        dZ(i) = 0.5 * dZ_2d(i,1) ! reset dZ to surface depth
        rhoSurf_2d(i,j) = rhoSurf(i)
        deltaRhoAtK(i) = 0.
        MLD(i,j) = 0.
        if (id_N2>0) then
          subMLN2(i,j) = 0.0
          dH_N2(i) = 0.0 ; dZ_N2(i) = 0.0
          T_subML(i) = 0.0  ; S_subML(i) = 0.0 ; T_deeper(i) = 0.0 ; S_deeper(i) = 0.0
          N2_region_set(i) = (G%mask2dT(i,j)<0.5) ! Only need to work on ocean points.
        endif
      enddo
    elseif (pRef_MLD(is) == 0.0) then
      rhoSurf(:) = 0.0
      do i=is,ie ; dZ(i) = 0.5 * dZ_2d(i,1) ; enddo ! Depth of center of surface layer
      call calculate_density(tv%T(:,j,1), tv%S(:,j,1), pRef_MLD, rhoSurf, tv%eqn_of_state, EOSdom)
      do i=is,ie
        rhoSurf_2d(i,j) = rhoSurf(i)
        deltaRhoAtK(i) = 0.
        MLD(i,j) = 0.
        if (id_N2>0) then
          subMLN2(i,j) = 0.0
          H_subML(i) = h(i,j,1) ; dH_N2(i) = 0.0 ; dZ_N2(i) = 0.0
          T_subML(i) = 0.0  ; S_subML(i) = 0.0 ; T_deeper(i) = 0.0 ; S_deeper(i) = 0.0
          N2_region_set(i) = (G%mask2dT(i,j)<0.5) ! Only need to work on ocean points.
        endif
      enddo
    endif

    do k=2,nz
      do i=is,ie
        dZm1(i) = dZ(i) ! Depth of center of layer K-1
        dZ(i) = dZ(i) + 0.5 * ( dZ_2d(i,k) + dZ_2d(i,k-1) ) ! Depth of center of layer K
      enddo

      ! Prepare to calculate stratification, N2, immediately below the mixed layer by finding
      ! the cells that extend over at least dz_subML.
      if (id_N2>0) then
         do i=is,ie
          if (MLD(i,j) == 0.0) then  ! Still in the mixed layer.
            H_subML(i) = H_subML(i) + h(i,j,k)
          elseif (.not.N2_region_set(i)) then ! This block is below the mixed layer, but N2 has not been found yet.
            if (dZ_N2(i) == 0.0) then ! Record the temperature, salinity, pressure, immediately below the ML
              T_subML(i) = tv%T(i,j,k) ; S_subML(i) = tv%S(i,j,k)
              H_subML(i) = H_subML(i) + 0.5 * h(i,j,k) ! Start midway through this layer.
              dH_N2(i) = 0.5 * h(i,j,k)
              dZ_N2(i) = 0.5 * dz_2d(i,k)
            elseif (dZ_N2(i) + dZ_2d(i,k) < dZ_sub_ML) then
              dH_N2(i) = dH_N2(i) + h(i,j,k)
              dZ_N2(i) = dZ_N2(i) + dz_2d(i,k)
            else  ! This layer includes the base of the region where N2 is calculated.
              T_deeper(i) = tv%T(i,j,k) ; S_deeper(i) = tv%S(i,j,k)
              dH_N2(i) = dH_N2(i) + 0.5 * h(i,j,k)
              dZ_N2(i) = dZ_N2(i) + 0.5 * dz_2d(i,k)
              N2_region_set(i) = .true.
            endif
          endif
        enddo ! i-loop
      endif ! id_N2>0

      ! Mixed-layer depth, using sigma-0 (surface reference pressure)
      do i=is,ie ; deltaRhoAtKm1(i) = deltaRhoAtK(i) ; enddo ! Store value from previous iteration of K
      call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pRef_MLD, deltaRhoAtK, tv%eqn_of_state, EOSdom)
      do i = is, ie
        deltaRhoAtK(i) = deltaRhoAtK(i) - rhoSurf(i) ! Density difference between layer K and surface
        ddRho = deltaRhoAtK(i) - deltaRhoAtKm1(i)
        if ((MLD(i,j) == 0.) .and. (ddRho > 0.) .and. &
            (deltaRhoAtKm1(i) < densityDiff) .and. (deltaRhoAtK(i) >= densityDiff)) then
          aFac = ( densityDiff - deltaRhoAtKm1(i) ) / ddRho
          MLD(i,j) = (dZ(i) * aFac + dZm1(i) * (1. - aFac))
        endif
        if (id_SQ > 0) MLD2(i,j) = MLD(i,j)**2
      enddo ! i-loop
    enddo ! k-loop
    do i=is,ie
      if ((MLD(i,j) == 0.) .and. (deltaRhoAtK(i) < densityDiff)) MLD(i,j) = dZ(i) ! Mixing goes to the bottom
    enddo

    if (id_N2>0) then  ! Now actually calculate stratification, N2, below the mixed layer.
      do i=is,ie ; pRef_N2(i) = (GV%g_Earth * GV%H_to_RZ) * (H_subML(i) + 0.5*dH_N2(i)) ; enddo
      ! if ((.not.N2_region_set(i)) .and. (dZ_N2(i) > 0.5*dZ_sub_ML)) then
      !    ! Use whatever stratification we can, measured over whatever distance is available?
      !    T_deeper(i) = tv%T(i,j,nz) ; S_deeper(i) = tv%S(i,j,nz)
      !    N2_region_set(i) = .true.
      ! endif
      call calculate_density(T_subML, S_subML, pRef_N2, rho_subML, tv%eqn_of_state, EOSdom)
      call calculate_density(T_deeper, S_deeper, pRef_N2, rho_deeper, tv%eqn_of_state, EOSdom)
      do i=is,ie ; if ((G%mask2dT(i,j) > 0.0) .and. N2_region_set(i)) then
        subMLN2(i,j) =  gE_rho0 * (rho_deeper(i) - rho_subML(i)) / dH_N2(i)
      endif ; enddo
    endif
  enddo ! j-loop

  if (id_MLD > 0) call post_data(id_MLD, MLD, diagPtr)
  if (id_N2 > 0)  call post_data(id_N2, subMLN2, diagPtr)
  if (id_SQ > 0)  call post_data(id_SQ, MLD2, diagPtr)

  if ((id_ref_z > 0) .and. (pRef_MLD(is)/=0.)) call post_data(id_ref_z, z_ref_diag , diagPtr)
  if (id_ref_rho > 0) call post_data(id_ref_rho, rhoSurf_2d , diagPtr)

end subroutine diagnoseMLDbyDensityDifference

!> Diagnose a mixed layer depth (MLD) determined by the depth a given energy value would mix.
!> This routine is appropriate in MOM_diabatic_aux due to its position within the time stepping.
subroutine diagnoseMLDbyEnergy(id_MLD, h, tv, G, GV, US, Mixing_Energy, diagPtr)
  ! Author: Brandon Reichl
  ! Date: October 2, 2020
  ! //
  ! *Note that gravity is assumed constant everywhere and divided out of all calculations.
  !
  ! This code has been written to step through the columns layer by layer, summing the PE
  ! change inferred by mixing the layer with all layers above.  When the change exceeds a
  ! threshold (determined by input array Mixing_Energy), the code needs to solve for how far
  ! into this layer the threshold PE change occurs (assuming constant density layers).
  ! This is expressed here via solving the function F(X) = 0 where:
  ! F(X) = 0.5 * ( Ca*X^3/(D1+X) + Cb*X^2/(D1+X) + Cc*X/(D1+X) + Dc/(D1+X)
  !                + Ca2*X^2 + Cb2*X + Cc2)
  ! where all coefficients are determined by the previous mixed layer depth, the
  ! density of the previous mixed layer, the present layer thickness, and the present
  ! layer density.  This equation is worked out by computing the total PE assuming constant
  ! density in the mixed layer as well as in the remaining part of the present layer that is
  ! not mixed.
  ! To solve for X in this equation a Newton's method iteration is employed, which
  ! converges extremely quickly (usually 1 guess) since this equation turns out to be rather
  ! linear for PE change with increasing X.
  ! Input parameters:
  integer, dimension(3),   intent(in) :: id_MLD      !< Energy output diagnostic IDs
  type(ocean_grid_type),   intent(in) :: G           !< Grid type
  type(verticalGrid_type), intent(in) :: GV          !< ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US          !< A dimensional unit scaling type
  real, dimension(3),      intent(in) :: Mixing_Energy !< Energy values for up to 3 MLDs [R Z3 T-2 ~> J m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h           !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in) :: tv          !< Structure containing pointers to any
                                                     !! available thermodynamic fields.
  type(diag_ctrl),         pointer    :: diagPtr     !< Diagnostics structure

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),3) :: MLD  ! Diagnosed mixed layer depth [Z ~> m].
  real, dimension(SZK_(GV)+1) :: Z_int    ! Depths of the interfaces from the surface [Z ~> m]
  real, dimension(SZI_(G),SZK_(GV)) :: dZ ! Layer thicknesses in depth units [Z ~> m]
  real, dimension(SZI_(G),SZK_(GV)) :: Rho_c ! Columns of layer densities [R ~> kg m-3]
  real, dimension(SZI_(G)) :: pRef_MLD   ! The reference pressure for the mixed layer
                                          ! depth calculation [R L2 T-2 ~> Pa]
  real, dimension(3) :: PE_threshold      ! The energy threshold divided by g [R Z2 ~> kg m-1]

  real :: PE_Threshold_fraction   ! The fractional tolerance of the specified energy
                                  ! for the energy used to mix to the diagnosed depth [nondim]
  real :: H_ML                    ! The accumulated depth of the mixed layer [Z ~> m]
  real :: PE                      ! The cumulative potential energy of the unmixed water column to a depth
                                  ! of H_ML, divided by the gravitational acceleration [R Z2 ~> kg m-1]
  real :: PE_Mixed                ! The potential energy of the completely mixed water column to a depth
                                  ! of H_ML, divided by the gravitational acceleration [R Z2 ~> kg m-1]
  real :: RhoDZ_ML                ! The depth integrated density of the mixed layer [R Z ~> kg m-2]
  real :: H_ML_TST                ! A new test value for the depth of the mixed layer [Z ~> m]
  real :: PE_Mixed_TST            ! The potential energy of the completely mixed water column to a depth
                                  ! of H_ML_TST, divided by the gravitational acceleration [R Z2 ~> kg m-1]
  real :: RhoDZ_ML_TST            ! A test value of the new depth integrated density of the mixed layer [R Z ~> kg m-2]
  real :: Rho_ML                  ! The average density of the mixed layer [R ~> kg m-3]

  ! These are all temporary variables used to shorten the expressions in the iterations.
  real :: R1, R2, Ca, Ca2 ! Some densities [R ~> kg m-3]
  real :: D1, D2, X, X2   ! Some thicknesses [Z ~> m]
  real :: Cb, Cb2    ! A depth integrated density [R Z ~> kg m-2]
  real :: C, D       ! A depth squared [Z2 ~> m2]
  real :: Cc, Cc2    ! A density times a depth squared [R Z2 ~> kg m-1]
  real :: Cd         ! A density times a depth cubed [R Z3 ~> kg]
  real :: Gx         ! A triple integral in depth of density [R Z3 ~> kg]
  real :: Gpx        ! The derivative of Gx with x  [R Z2 ~> kg m-1]
  real :: Hx         ! The vertical integral depth [Z ~> m]
  real :: iHx        ! The inverse of Hx [Z-1 ~> m-1]
  real :: Hpx        ! The derivative of Hx with x, hard coded to 1.  Why? [nondim]
  real :: Ix         ! A double integral in depth of density [R Z2 ~> kg m-1]
  real :: Ipx        ! The derivative of Ix with x [R Z ~> kg m-2]
  real :: Fgx        ! The mixing energy difference from the target [R Z2 ~> kg m-1]
  real :: Fpx        ! The derivative of Fgx with x  [R Z ~> kg m-2]

  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: IT, iM
  integer :: i, j, is, ie, js, je, k, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  pRef_MLD(:) = 0.0
  mld(:,:,:) = 0.0
  PE_Threshold_fraction = 1.e-4 !Fixed threshold of 0.01%, could be runtime.

  do iM=1,3
    PE_threshold(iM) = Mixing_Energy(iM) / (US%L_to_Z**2*GV%g_Earth)
  enddo

  MLD(:,:,:) = 0.0

  EOSdom(:) = EOS_domain(G%HI)

  do j=js,je
    ! Find the vertical distances across layers.
    call thickness_to_dz(h, tv, dz, j, G, GV)

    do k=1,nz
      call calculate_density(tv%T(:,j,k), tv%S(:,j,K), pRef_MLD, rho_c(:,k), tv%eqn_of_state, EOSdom)
    enddo

    do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then

      Z_int(1) = 0.0
      do k=1,nz
        Z_int(K+1) = Z_int(K) - dZ(i,k)
      enddo

      do iM=1,3

        ! Initialize these for each column-wise calculation
        PE = 0.0
        RhoDZ_ML = 0.0
        H_ML = 0.0
        RhoDZ_ML_TST = 0.0
        H_ML_TST = 0.0
        PE_Mixed = 0.0

        do k=1,nz

          ! This is the unmixed PE cumulative sum from top down
          PE = PE + 0.5 * Rho_c(i,k) * (Z_int(K)**2 - Z_int(K+1)**2)

          ! This is the depth and integral of density
          H_ML_TST = H_ML + dZ(i,k)
          RhoDZ_ML_TST = RhoDZ_ML + Rho_c(i,k) * dZ(i,k)

          ! The average density assuming all layers including this were mixed
          Rho_ML = RhoDZ_ML_TST/H_ML_TST

          ! The PE assuming all layers including this were mixed
          ! Note that 0. could be replaced with "Surface", which doesn't have to be 0
          ! but 0 is a good reference value.
          PE_Mixed_TST = 0.5 * Rho_ML * (0.**2 - (0. - H_ML_TST)**2)

          ! Check if we supplied enough energy to mix to this layer
          if (PE_Mixed_TST - PE <= PE_threshold(iM)) then
            H_ML = H_ML_TST
            RhoDZ_ML = RhoDZ_ML_TST

          else ! If not, we need to solve where the energy ran out
            ! This will be done with a Newton's method iteration:

            R1 = RhoDZ_ML / H_ML ! The density of the mixed layer (not including this layer)
            D1 = H_ML ! The thickness of the mixed layer (not including this layer)
            R2 = Rho_c(i,k) ! The density of this layer
            D2 = dZ(i,k) ! The thickness of this layer

            ! This block could be used to calculate the function coefficients if
            ! we don't reference all values to a surface designated as z=0
            ! S = Surface
            ! Ca  = -(R2)
            ! Cb  = -( (R1*D1) + R2*(2.*D1-2.*S) )
            ! D   = D1**2. - 2.*D1*S
            ! Cc  = -( R1*D1*(2.*D1-2.*S) + R2*D )
            ! Cd  = -(R1*D1*D)
            ! Ca2 = R2
            ! Cb2 = R2*(2*D1-2*S)
            ! C   = S**2 + D2**2 + D1**2 - 2*D1*S - 2.*D2*S +2.*D1*D2
            ! Cc2 = R2*(D+S**2-C)
            !
            ! If the surface is S = 0, it simplifies to:
            Ca  = -R2
            Cb  = -(R1 * D1 + R2 * (2. * D1))
            D   = D1**2
            Cc  = -(R1 * D1 * (2. * D1) + (R2 * D))
            Cd  = -R1 * (D1 * D)
            Ca2 = R2
            Cb2 = R2 * (2. * D1)
            C   = D2**2 + D1**2 + 2. * (D1 * D2)
            Cc2 = R2 * (D - C)

            ! First guess for an iteration using Newton's method
            X = dZ(i,k) * 0.5

            IT=0
            do while(IT<10)!We can iterate up to 10 times
              ! We are trying to solve the function:
              ! F(x) = G(x)/H(x)+I(x)
              ! for where F(x) = PE+PE_threshold, or equivalently for where
              ! F(x) = G(x)/H(x)+I(x) - (PE+PE_threshold) = 0
              ! We also need the derivative of this function for the Newton's method iteration
              ! F'(x) = (G'(x)H(x)-G(x)H'(x))/H(x)^2 + I'(x)
              ! G and its derivative
              Gx = 0.5 * (Ca * (X*X*X) + Cb * X**2 + Cc * X + Cd)
              Gpx = 0.5 * (3. * (Ca * X**2) + 2. * (Cb * X) + Cc)
              ! H, its inverse, and its derivative
              Hx = D1 + X
              iHx = 1. / Hx
              Hpx = 1.
              ! I and its derivative
              Ix = 0.5 * (Ca2 * X**2 + Cb2 * X + Cc2)
              Ipx = 0.5 * (2. * Ca2 * X + Cb2)

              ! The Function and its derivative:
              PE_Mixed = Gx * iHx + Ix
              Fgx = PE_Mixed - (PE + PE_threshold(iM))
              Fpx = (Gpx * Hx - Hpx * Gx) * iHx**2 + Ipx

              ! Check if our solution is within the threshold bounds, if not update
              ! using Newton's method.  This appears to converge almost always in
              ! one step because the function is very close to linear in most applications.
              if (abs(Fgx) > PE_Threshold(iM) * PE_Threshold_fraction) then
                X2 = X - Fgx / Fpx
                IT = IT + 1
                if (X2 < 0. .or. X2 > dZ(i,k)) then
                  ! The iteration seems to be robust, but we need to do something *if*
                  ! things go wrong... How should we treat failed iteration?
                  ! Present solution: Stop trying to compute and just say we can't mix this layer.
                  X=0
                  exit
                else
                  X = X2
                endif
              else
                exit! Quit the iteration
              endif
            enddo
            H_ML = H_ML + X
            exit! Quit looping through the column
          endif
        enddo
        MLD(i,j,iM) = H_ML
      enddo
    endif ; enddo
  enddo

  if (id_MLD(1) > 0) call post_data(id_MLD(1), MLD(:,:,1), diagPtr)
  if (id_MLD(2) > 0) call post_data(id_MLD(2), MLD(:,:,2), diagPtr)
  if (id_MLD(3) > 0) call post_data(id_MLD(3), MLD(:,:,3), diagPtr)

end subroutine diagnoseMLDbyEnergy

!> \namespace mom_diagnose_mld
!!
!!    This module contains subroutines that apply various diabatic processes.  Usually these
!!  subroutines are called from the MOM_diabatic module.  All of these routines use appropriate
!!  limiters or logic to work properly with arbitrary layer thicknesses (including massless layers)
!!  and an arbitrarily large timestep.
!!
!!    The subroutine diagnoseMLDbyDensityDifference diagnoses a mixed layer depth based on a
!!  density difference criterion, and may also estimate the stratification of the water below
!!  this diagnosed mixed layer.
!!
!!    The subroutine diagnoseMLDbyEnergy diagnoses a mixed layer depth based on a mixing-energy
!!  criterion, as described by Reichl et al., 2022, JGR: Oceans, doi:10.1029/2021JC018140.
!!

end module MOM_diagnose_mld
