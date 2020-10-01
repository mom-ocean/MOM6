!> Implements the thermodynamic aspects of ocean / ice-shelf interactions,
!!  along with a crude placeholder for a later implementation of full
!!  ice shelf dynamics, all using the MOM framework and coding style.
module MOM_ice_shelf_state

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_ROUTINE
use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_grid, only : MOM_grid_init, ocean_grid_type
use MOM_get_input, only : directories, Get_MOM_input
use mpp_mod, only : mpp_sum, mpp_max, mpp_min, mpp_pe, mpp_npes, mpp_sync
use MOM_coms, only : reproducing_sum
use MOM_checksums, only : hchksum, qchksum, chksum, uchksum, vchksum, uvchksum

implicit none ; private

public ice_shelf_state_end, ice_shelf_state_init

!> Structure that describes the ice shelf state
type, public :: ice_shelf_state
  real, pointer, dimension(:,:) :: &
    mass_shelf => NULL(), &    !< The mass per unit area of the ice shelf or sheet [R Z ~> kg m-2].
    area_shelf_h => NULL(), &  !< The area per cell covered by the ice shelf [L2 ~> m2].
    h_shelf => NULL(), &       !< the thickness of the shelf [Z ~> m], redundant with mass but may
                               !! make the code more readable
    hmask => NULL(),&          !< Mask used to indicate ice-covered or partiall-covered cells
                               !! 1: fully covered, solve for velocity here (for now all
                               !!   ice-covered cells are treated the same, this may change)
                               !! 2: partially covered, do not solve for velocity
                               !! 0: no ice in cell.
                               !! 3: bdry condition on thickness set - not in computational domain
                               !! -2 : default (out of computational boundary, and) not = 3
                               !! NOTE: hmask will change over time and NEEDS TO BE MAINTAINED
                               !!   otherwise the wrong nodes will be included in velocity calcs.

    tflux_ocn => NULL(), &     !< The downward sensible ocean heat flux at the
                               !! ocean-ice interface [Q R Z T-1 ~> W m-2].
    salt_flux => NULL(), &     !< The downward salt flux at the ocean-ice
                               !! interface [kgSalt kgWater-1 R Z T-1 ~> kgSalt m-2 s-1].
    water_flux => NULL(), &    !< The net downward liquid water flux at the
                               !! ocean-ice interface [R Z T-1 ~> kg m-2 s-1].
    tflux_shelf => NULL(), &   !< The downward diffusive heat flux in the ice
                               !! shelf at the ice-ocean interface [Q R Z T-1 ~> W m-2].

    tfreeze => NULL()          !< The freezing point potential temperature
                               !! an the ice-ocean interface [degC].

end type ice_shelf_state

contains

!> Deallocates all memory associated with this module
subroutine ice_shelf_state_init(ISS, G)
  type(ice_shelf_state), pointer    :: ISS !< A pointer to the ice shelf state structure
  type(ocean_grid_type), intent(in) :: G   !< The grid structure used by the ice shelf.

  integer :: isd, ied, jsd, jed
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed

  if (associated(ISS)) then
    call MOM_error(FATAL, "MOM_ice_shelf_state.F90, ice_shelf_state_init: "// &
                          "called with an associated ice_shelf_state pointer.")
    return
  endif
  allocate(ISS)

  allocate(ISS%mass_shelf(isd:ied,jsd:jed) )   ; ISS%mass_shelf(:,:) = 0.0
  allocate(ISS%area_shelf_h(isd:ied,jsd:jed) ) ; ISS%area_shelf_h(:,:) = 0.0
  allocate(ISS%h_shelf(isd:ied,jsd:jed) )      ; ISS%h_shelf(:,:) = 0.0
  allocate(ISS%hmask(isd:ied,jsd:jed) )        ; ISS%hmask(:,:) = -2.0

  allocate(ISS%tflux_ocn(isd:ied,jsd:jed) )    ; ISS%tflux_ocn(:,:) = 0.0
  allocate(ISS%water_flux(isd:ied,jsd:jed) )   ; ISS%water_flux(:,:) = 0.0
  allocate(ISS%salt_flux(isd:ied,jsd:jed) )    ; ISS%salt_flux(:,:) = 0.0
  allocate(ISS%tflux_shelf(isd:ied,jsd:jed) )  ; ISS%tflux_shelf(:,:) = 0.0
  allocate(ISS%tfreeze(isd:ied,jsd:jed) )      ; ISS%tfreeze(:,:) = 0.0

end subroutine ice_shelf_state_init


!> Deallocates all memory associated with this module
subroutine ice_shelf_state_end(ISS)
  type(ice_shelf_state), pointer :: ISS !< A pointer to the ice shelf state structure

  if (.not.associated(ISS)) return

  deallocate(ISS%mass_shelf, ISS%area_shelf_h, ISS%h_shelf, ISS%hmask)

  deallocate(ISS%tflux_ocn, ISS%water_flux, ISS%salt_flux, ISS%tflux_shelf)
  deallocate(ISS%tfreeze)

  deallocate(ISS)

end subroutine ice_shelf_state_end


end module MOM_ice_shelf_state
