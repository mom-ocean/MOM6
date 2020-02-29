!> Routines incorporating the effects of marine ice (sea-ice and icebergs) into
!! the ocean model dynamics and thermodynamics.
module MOM_marine_ice

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_constants,     only : hlf
use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_domains,       only : TO_ALL, Omit_Corners
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_forcing_type,  only : allocate_forcing_type
use MOM_forcing_type,  only : forcing, mech_forcing
use MOM_grid,          only : ocean_grid_type
use MOM_time_manager,  only : time_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : surface

implicit none ; private

#include <MOM_memory.h>

public iceberg_forces, iceberg_fluxes, marine_ice_init

!> Control structure for MOM_marine_ice
type, public :: marine_ice_CS ; private
  real :: kv_iceberg          !< The viscosity of the icebergs [m2 s-1] (for ice rigidity)
  real :: berg_area_threshold !< Fraction of grid cell which iceberg must occupy
                              !! so that fluxes below are set to zero. (0.5 is a
                              !! good value to use.) Not applied for negative values.
  real :: latent_heat_fusion  !< Latent heat of fusion [J kg-1]
  real :: density_iceberg     !< A typical density of icebergs [kg m-3] (for ice rigidity)

  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the timing of diagnostic output.
end type marine_ice_CS

contains

!> add_berg_flux_to_shelf adds rigidity and ice-area coverage due to icebergs
!! to the forces type fields, and adds ice-areal coverage and modifies various
!! thermodynamic fluxes due to the presence of icebergs.
subroutine iceberg_forces(G, forces, use_ice_shelf, sfc_state, &
                                  time_step, CS)
  type(ocean_grid_type), intent(inout) :: G       !< The ocean's grid structure
  type(mech_forcing),    intent(inout) :: forces  !< A structure with the driving mechanical forces
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  logical,               intent(in)    :: use_ice_shelf  !< If true, this configuration uses ice shelves.
  real,                  intent(in)    :: time_step   !< The coupling time step [s].
  type(marine_ice_CS),   pointer       :: CS !< Pointer to the control structure for MOM_marine_ice

  real :: kv_rho_ice ! The viscosity of ice divided by its density [m5 kg-1 s-1].
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  !This routine adds iceberg data to the ice shelf data (if ice shelf is used)
  !which can then be used to change the top of ocean boundary condition used in
  !the ocean model. This routine is taken from the add_shelf_flux subroutine
  !within the ice shelf model.

  if (.not.associated(CS)) return

  if (.not.(associated(forces%area_berg) .and.  associated(forces%mass_berg) ) ) return

  if (.not.(associated(forces%frac_shelf_u) .and. associated(forces%frac_shelf_v) .and. &
            associated(forces%rigidity_ice_u) .and. associated(forces%rigidity_ice_v)) ) return

  ! This section sets or augments the values of fields in forces.
  if (.not. use_ice_shelf) then
    forces%frac_shelf_u(:,:) = 0.0 ; forces%frac_shelf_v(:,:) = 0.0
  endif
  if (.not. forces%accumulate_rigidity) then
    forces%rigidity_ice_u(:,:) = 0.0 ; forces%rigidity_ice_v(:,:) = 0.0
  endif

  call pass_var(forces%area_berg, G%domain, TO_ALL+Omit_corners, halo=1, complete=.false.)
  call pass_var(forces%mass_berg, G%domain, TO_ALL+Omit_corners, halo=1, complete=.true.)
  kv_rho_ice = CS%kv_iceberg / CS%density_iceberg
  do j=js,je ; do I=is-1,ie
    if ((G%areaT(i,j) + G%areaT(i+1,j) > 0.0)) & ! .and. (G%dxdy_u(I,j) > 0.0)) &
      forces%frac_shelf_u(I,j) = forces%frac_shelf_u(I,j) + &
          (((forces%area_berg(i,j)*G%US%L_to_m**2*G%areaT(i,j)) + &
            (forces%area_berg(i+1,j)*G%US%L_to_m**2*G%areaT(i+1,j))) / &
           (G%US%L_to_m**2*G%areaT(i,j) + G%US%L_to_m**2*G%areaT(i+1,j)) )
    forces%rigidity_ice_u(I,j) = forces%rigidity_ice_u(I,j) + kv_rho_ice * &
                         min(forces%mass_berg(i,j), forces%mass_berg(i+1,j))
  enddo ; enddo
  do J=js-1,je ; do i=is,ie
    if ((G%areaT(i,j) + G%areaT(i,j+1) > 0.0)) & ! .and. (G%dxdy_v(i,J) > 0.0)) &
      forces%frac_shelf_v(i,J) = forces%frac_shelf_v(i,J) + &
          (((forces%area_berg(i,j)*G%US%L_to_m**2*G%areaT(i,j)) + &
            (forces%area_berg(i,j+1)*G%US%L_to_m**2*G%areaT(i,j+1))) / &
           (G%US%L_to_m**2*G%areaT(i,j) + G%US%L_to_m**2*G%areaT(i,j+1)) )
    forces%rigidity_ice_v(i,J) = forces%rigidity_ice_v(i,J) + kv_rho_ice * &
                         min(forces%mass_berg(i,j), forces%mass_berg(i,j+1))
  enddo ; enddo
  !### This halo update may be unnecessary. Test it.  -RWH
  call pass_vector(forces%frac_shelf_u, forces%frac_shelf_v, G%domain, TO_ALL, CGRID_NE)

end subroutine iceberg_forces

!> iceberg_fluxes adds ice-area-coverage and modifies various
!! thermodynamic fluxes due to the presence of icebergs.
subroutine iceberg_fluxes(G, US, fluxes, use_ice_shelf, sfc_state, &
                          time_step, CS)
  type(ocean_grid_type), intent(inout) :: G       !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US      !< A dimensional unit scaling type
  type(forcing),         intent(inout) :: fluxes  !< A structure with pointers to themodynamic,
                                                  !! tracer and mass exchange forcing fields
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  logical,               intent(in)    :: use_ice_shelf  !< If true, this configuration uses ice shelves.
  real,                  intent(in)    :: time_step   !< The coupling time step [s].
  type(marine_ice_CS),   pointer       :: CS !< Pointer to the control structure for MOM_marine_ice

  real :: fraz      ! refreezing rate [R Z T-1 ~> kg m-2 s-1]
  real :: I_dt_LHF  ! The inverse of the timestep times the latent heat of fusion [kg J-1 T-1 ~> kg J-1 s-1].
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
  !This routine adds iceberg data to the ice shelf data (if ice shelf is used)
  !which can then be used to change the top of ocean boundary condition used in
  !the ocean model. This routine is taken from the add_shelf_flux subroutine
  !within the ice shelf model.

  if (.not.associated(CS)) return
  if (.not.(associated(fluxes%area_berg) .and. associated(fluxes%ustar_berg) .and. &
            associated(fluxes%mass_berg) ) ) return
  if (.not.(associated(fluxes%frac_shelf_h) .and. associated(fluxes%ustar_shelf)) ) return


  if (.not.(associated(fluxes%area_berg) .and. associated(fluxes%ustar_berg) .and. &
            associated(fluxes%mass_berg) ) ) return
  if (.not. use_ice_shelf) then
    fluxes%frac_shelf_h(:,:) = 0.
    fluxes%ustar_shelf(:,:) = 0.
  endif
  do j=jsd,jed ; do i=isd,ied ; if (G%areaT(i,j) > 0.0) then
    fluxes%frac_shelf_h(i,j) = fluxes%frac_shelf_h(i,j) + fluxes%area_berg(i,j)
    fluxes%ustar_shelf(i,j)  = fluxes%ustar_shelf(i,j)  + fluxes%ustar_berg(i,j)
  endif ; enddo ; enddo

  !Zero'ing out other fluxes under the tabular icebergs
  if (CS%berg_area_threshold >= 0.) then
    I_dt_LHF = 1.0 / (US%s_to_T*time_step * CS%latent_heat_fusion)
    do j=jsd,jed ; do i=isd,ied
      if (fluxes%frac_shelf_h(i,j) > CS%berg_area_threshold) then
        ! Only applying for ice shelf covering most of cell.

        if (associated(fluxes%sw)) fluxes%sw(i,j) = 0.0
        if (associated(fluxes%lw)) fluxes%lw(i,j) = 0.0
        if (associated(fluxes%latent)) fluxes%latent(i,j) = 0.0
        if (associated(fluxes%evap)) fluxes%evap(i,j) = 0.0

        ! Add frazil formation diagnosed by the ocean model [J m-2] in the
        ! form of surface layer evaporation [R Z T-1 ~> kg m-2 s-1]. Update lprec in the
        ! control structure for diagnostic purposes.

        if (allocated(sfc_state%frazil)) then
          fraz = US%kg_m3_to_R*US%m_to_Z*sfc_state%frazil(i,j) * I_dt_LHF
          if (associated(fluxes%evap))  fluxes%evap(i,j)  = fluxes%evap(i,j)  - fraz
        ! if (associated(fluxes%lprec)) fluxes%lprec(i,j) = fluxes%lprec(i,j) - fraz
          sfc_state%frazil(i,j) = 0.0
        endif

        !Alon: Should these be set to zero too?
        if (associated(fluxes%sens)) fluxes%sens(i,j) = 0.0
        if (associated(fluxes%salt_flux)) fluxes%salt_flux(i,j) = 0.0
        if (associated(fluxes%lprec)) fluxes%lprec(i,j) = 0.0
      endif
    enddo ; enddo
  endif

end subroutine iceberg_fluxes

!> Initialize control structure for MOM_marine_ice
subroutine marine_ice_init(Time, G, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G !< Ocean grid structure
  type(param_file_type),   intent(in)    :: param_file !< Runtime parameter handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(marine_ice_CS),     pointer       :: CS   !< Pointer to the control structure for MOM_marine_ice
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_marine_ice"  ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "marine_ice_init called with an associated control structure.")
    return
  else ; allocate(CS) ; endif

  ! Write all relevant parameters to the model log.
  call log_version(mdl, version)

  call get_param(param_file, mdl, "KV_ICEBERG",  CS%kv_iceberg, &
                 "The viscosity of the icebergs",  units="m2 s-1",default=1.0e10)
  call get_param(param_file, mdl, "DENSITY_ICEBERGS",  CS%density_iceberg, &
                 "A typical density of icebergs.", units="kg m-3", default=917.0)
  call get_param(param_file, mdl, "LATENT_HEAT_FUSION", CS%latent_heat_fusion, &
                 "The latent heat of fusion.", units="J/kg", default=hlf)
  call get_param(param_file, mdl, "BERG_AREA_THRESHOLD", CS%berg_area_threshold, &
                 "Fraction of grid cell which iceberg must occupy, so that fluxes "//&
                 "below berg are set to zero. Not applied for negative "//&
                 "values.", units="non-dim", default=-1.0)

end subroutine marine_ice_init

end module MOM_marine_ice
