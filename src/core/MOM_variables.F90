!> Provides transparent structures with groups of MOM6 variables and supporting routines
module MOM_variables

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_array_transform, only : rotate_array, rotate_vector
use MOM_coupler_types, only : coupler_1d_bc_type, coupler_2d_bc_type
use MOM_coupler_types, only : coupler_type_spawn, coupler_type_destructor, coupler_type_initialized
use MOM_debugging,     only : hchksum
use MOM_domains,       only : MOM_domain_type, get_domain_extent, group_pass_type
use MOM_EOS,           only : EOS_type
use MOM_error_handler, only : MOM_error, FATAL
use MOM_grid,          only : ocean_grid_type
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public allocate_surface_state, deallocate_surface_state, MOM_thermovar_chksum
public ocean_grid_type, alloc_BT_cont_type, dealloc_BT_cont_type
public rotate_surface_state

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> A structure for creating arrays of pointers to 3D arrays
type, public :: p3d
  real, dimension(:,:,:), pointer :: p => NULL() !< A pointer to a 3D array
end type p3d
!> A structure for creating arrays of pointers to 2D arrays
type, public :: p2d
  real, dimension(:,:), pointer :: p => NULL() !< A pointer to a 2D array
end type p2d

!> Pointers to various fields which may be used describe the surface state of MOM, and which
!! will be returned to a the calling program
type, public :: surface
  real, allocatable, dimension(:,:) :: &
    SST, &         !< The sea surface temperature [degC].
    SSS, &         !< The sea surface salinity [ppt ~> psu or gSalt/kg].
    sfc_density, & !< The mixed layer density [R ~> kg m-3].
    Hml, &         !< The mixed layer depth [Z ~> m].
    u, &           !< The mixed layer zonal velocity [L T-1 ~> m s-1].
    v, &           !< The mixed layer meridional velocity [L T-1 ~> m s-1].
    sea_lev, &     !< The sea level [Z ~> m].  If a reduced surface gravity is
                   !! used, that is compensated for in sea_lev.
    frazil, &      !< The energy needed to heat the ocean column to the freezing point during
                   !! the call to step_MOM [Q R Z ~> J m-2].
    melt_potential, & !< Instantaneous amount of heat that can be used to melt sea ice [Q R Z ~> J m-2].
                      !! This is computed w.r.t. surface freezing temperature.
    ocean_mass, &  !< The total mass of the ocean [R Z ~> kg m-2].
    ocean_heat, &  !< The total heat content of the ocean in [degC R Z ~> degC kg m-2].
    ocean_salt, &  !< The total salt content of the ocean in [kgSalt kg-1 R Z ~> kgSalt m-2].
    taux_shelf, &  !< The zonal stresses on the ocean under shelves [R L Z T-2 ~> Pa].
    tauy_shelf, &  !< The meridional stresses on the ocean under shelves [R L Z T-2 ~> Pa].
    TempxPmE, &    !< The net inflow of water into the ocean times the temperature at which this
                   !! inflow occurs during the call to step_MOM [degC R Z ~> degC kg m-2].
    salt_deficit, & !< The salt needed to maintain the ocean column above a minimum
                   !! salinity over the call to step_MOM [kgSalt kg-1 R Z ~> kgSalt m-2].
    internal_heat  !< Any internal or geothermal heat sources that are applied to the ocean
                   !! integrated over the call to step_MOM [degC R Z ~> degC kg m-2].
  logical :: T_is_conT = .false. !< If true, the temperature variable SST is actually the
                   !! conservative temperature in [degC].
  logical :: S_is_absS = .false. !< If true, the salinity variable SSS is actually the
                   !! absolute salinity in [g/kg].
  type(coupler_2d_bc_type) :: tr_fields !< A structure that may contain an
                !! array of named fields describing tracer-related quantities.
       !### NOTE: ALL OF THE ARRAYS IN TR_FIELDS USE THE COUPLER'S INDEXING CONVENTION AND HAVE NO
       !###       HALOS!  THIS IS DONE TO CONFORM TO THE TREATMENT IN MOM4, BUT I DON'T LIKE IT! -RWH
  logical :: arrays_allocated = .false.  !< A flag that indicates whether the surface type
                !! has had its memory allocated.
end type surface

!> Pointers to an assortment of thermodynamic fields that may be available, including
!! potential temperature, salinity, heat capacity, and the equation of state control structure.
type, public :: thermo_var_ptrs
  ! If allocated, the following variables have nz layers.
  real, pointer :: T(:,:,:) => NULL() !< Potential temperature [degC].
  real, pointer :: S(:,:,:) => NULL() !< Salinity [PSU] or [gSalt/kg], generically [ppt].
  real, pointer :: p_surf(:,:) => NULL() !< Ocean surface pressure used in equation of state
                         !! calculations [R L2 T-2 ~> Pa]
  type(EOS_type), pointer :: eqn_of_state => NULL() !< Type that indicates the
                         !! equation of state to use.
  real :: P_Ref          !<   The coordinate-density reference pressure [R L2 T-2 ~> Pa].
                         !! This is the pressure used to calculate Rml from
                         !! T and S when eqn_of_state is associated.
  real :: C_p            !<   The heat capacity of seawater [Q degC-1 ~> J degC-1 kg-1].
                         !! When conservative temperature is used, this is
                         !! constant and exactly 3991.86795711963 J degC-1 kg-1.
  logical :: T_is_conT = .false. !< If true, the temperature variable tv%T is
                         !! actually the conservative temperature [degC].
  logical :: S_is_absS = .false. !< If true, the salinity variable tv%S is
                         !! actually the absolute salinity in units of [gSalt/kg].
  real :: min_salinity = 0.01 !< The minimum value of salinity when BOUND_SALINITY=True [ppt].
                         !! The default is 0.01 for backward compatibility but should be 0.
  ! These arrays are accumulated fluxes for communication with other components.
  real, dimension(:,:), pointer :: frazil => NULL()
                         !< The energy needed to heat the ocean column to the
                         !! freezing point since calculate_surface_state was2
                         !! last called [Q Z R ~> J m-2].
  real, dimension(:,:), pointer :: salt_deficit => NULL()
                         !<   The salt needed to maintain the ocean column
                         !! at a minimum salinity of MIN_SALINITY since the last time
                         !! that calculate_surface_state was called, [ppt R Z ~> gSalt m-2].
  real, dimension(:,:), pointer :: TempxPmE => NULL()
                         !<   The net inflow of water into the ocean times the
                         !! temperature at which this inflow occurs since the
                         !! last call to calculate_surface_state [degC R Z ~> degC kg m-2].
                         !! This should be prescribed in the forcing fields, but
                         !! as it often is not, this is a useful heat budget diagnostic.
  real, dimension(:,:), pointer :: internal_heat => NULL()
                         !< Any internal or geothermal heat sources that
                         !! have been applied to the ocean since the last call to
                         !! calculate_surface_state [degC R Z ~> degC kg m-2].
  ! The following variables are most normally not used but when they are they
  ! will be either set by parameterizations or prognostic.
  real, pointer :: varT(:,:,:) => NULL() !< SGS variance of potential temperature [degC2].
  real, pointer :: varS(:,:,:) => NULL() !< SGS variance of salinity [ppt2].
  real, pointer :: covarTS(:,:,:) => NULL() !< SGS covariance of salinity and potential
                                  !! temperature [degC ppt].
end type thermo_var_ptrs

!> Pointers to all of the prognostic variables allocated in MOM_variables.F90 and MOM.F90.
!!
!! It is useful for sending these variables for diagnostics, and in preparation for ensembles
!! later on.  All variables have the same names as the local (public) variables
!! they refer to in MOM.F90.
type, public :: ocean_internal_state
  real, pointer, dimension(:,:,:) :: &
    T => NULL(), & !< Pointer to the temperature state variable [degC]
    S => NULL(), & !< Pointer to the salinity state variable [ppt ~> PSU or g/kg]
    u => NULL(), & !< Pointer to the zonal velocity [L T-1 ~> m s-1]
    v => NULL(), & !< Pointer to the meridional velocity [L T-1 ~> m s-1]
    h => NULL()    !< Pointer to the layer thicknesses [H ~> m or kg m-2]
  real, pointer, dimension(:,:,:) :: &
    uh => NULL(), & !<  Pointer to zonal transports [H L2 T-1 ~> m3 s-1 or kg s-1]
    vh => NULL()    !<  Pointer to meridional transports [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, pointer, dimension(:,:,:) :: &
    CAu => NULL(), & !< Pointer to the zonal Coriolis and Advective acceleration [L T-2 ~> m s-2]
    CAv => NULL(), & !< Pointer to the meridional Coriolis and Advective acceleration [L T-2 ~> m s-2]
    PFu => NULL(), & !< Pointer to the zonal Pressure force acceleration [L T-2 ~> m s-2]
    PFv => NULL(), & !< Pointer to the meridional Pressure force acceleration [L T-2 ~> m s-2]
    diffu => NULL(), & !< Pointer to the zonal acceleration due to lateral viscosity [L T-2 ~> m s-2]
    diffv => NULL(), & !< Pointer to the meridional acceleration due to lateral viscosity [L T-2 ~> m s-2]
    pbce => NULL(), &  !< Pointer to the baroclinic pressure force dependency on free surface movement
                       !! [L2 T-2 H-1 ~> m s-2 or m4 kg-1 s-2]
    u_accel_bt => NULL(), & !< Pointer to the zonal barotropic-solver acceleration [L T-2 ~> m s-2]
    v_accel_bt => NULL()  !< Pointer to the meridional barotropic-solver acceleration [L T-2 ~> m s-2]
  real, pointer, dimension(:,:,:) :: &
    u_av => NULL(), &  !< Pointer to zonal velocity averaged over the timestep [L T-1 ~> m s-1]
    v_av => NULL(), &  !< Pointer to meridional velocity averaged over the timestep [L T-1 ~> m s-1]
    u_prev => NULL(), & !< Pointer to zonal velocity at the end of the last timestep [L T-1 ~> m s-1]
    v_prev => NULL()   !< Pointer to meridional velocity at the end of the last timestep [L T-1 ~> m s-1]
end type ocean_internal_state

!> Pointers to arrays with accelerations, which can later be used for derived diagnostics, like energy balances.
type, public :: accel_diag_ptrs

  ! Each of the following fields has nz layers.
  real, pointer, dimension(:,:,:) :: &
    diffu => NULL(), &     !< Zonal acceleration due to along isopycnal viscosity [L T-2 ~> m s-2]
    diffv => NULL(), &     !< Meridional acceleration due to along isopycnal viscosity [L T-2 ~> m s-2]
    CAu => NULL(), &       !< Zonal Coriolis and momentum advection accelerations [L T-2 ~> m s-2]
    CAv => NULL(), &       !< Meridional Coriolis and momentum advection accelerations [L T-2 ~> m s-2]
    PFu => NULL(), &       !< Zonal acceleration due to pressure forces [L T-2 ~> m s-2]
    PFv => NULL(), &       !< Meridional acceleration due to pressure forces [L T-2 ~> m s-2]
    du_dt_visc => NULL(), &!< Zonal acceleration due to vertical viscosity [L T-2 ~> m s-2]
    dv_dt_visc => NULL(), &!< Meridional acceleration due to vertical viscosity [L T-2 ~> m s-2]
    du_dt_dia => NULL(), & !< Zonal acceleration due to diapycnal  mixing [L T-2 ~> m s-2]
    dv_dt_dia => NULL(), & !< Meridional acceleration due to diapycnal  mixing [L T-2 ~> m s-2]
    u_accel_bt => NULL(), &!< Pointer to the zonal barotropic-solver acceleration [L T-2 ~> m s-2]
    v_accel_bt => NULL()   !< Pointer to the meridional barotropic-solver acceleration [L T-2 ~> m s-2]
  real, pointer, dimension(:,:,:) :: du_other => NULL()
                           !< Zonal velocity changes due to any other processes that are
                           !! not due to any explicit accelerations [L T-1 ~> m s-1].
  real, pointer, dimension(:,:,:) :: dv_other => NULL()
                           !< Meridional velocity changes due to any other processes that are
                           !! not due to any explicit accelerations [L T-1 ~> m s-1].

  ! These accelerations are sub-terms included in the accelerations above.
  real, pointer :: gradKEu(:,:,:) => NULL()  !< gradKEu = - d/dx(u2) [L T-2 ~> m s-2]
  real, pointer :: gradKEv(:,:,:) => NULL()  !< gradKEv = - d/dy(u2) [L T-2 ~> m s-2]
  real, pointer :: rv_x_v(:,:,:) => NULL()   !< rv_x_v = rv * v at u [L T-2 ~> m s-2]
  real, pointer :: rv_x_u(:,:,:) => NULL()   !< rv_x_u = rv * u at v [L T-2 ~> m s-2]

  real, pointer :: diag_hfrac_u(:,:,:) => NULL() !< Fractional layer thickness at u points
  real, pointer :: diag_hfrac_v(:,:,:) => NULL() !< Fractional layer thickness at v points
  real, pointer :: diag_hu(:,:,:) => NULL() !< layer thickness at u points
  real, pointer :: diag_hv(:,:,:) => NULL() !< layer thickness at v points

end type accel_diag_ptrs

!> Pointers to arrays with transports, which can later be used for derived diagnostics, like energy balances.
type, public :: cont_diag_ptrs

! Each of the following fields has nz layers.
  real, pointer, dimension(:,:,:) :: &
    uh => NULL(), &   !< Resolved zonal layer thickness fluxes, [H L2 T-1 ~> m3 s-1 or kg s-1]
    vh => NULL(), &   !< Resolved meridional layer thickness fluxes, [H L2 T-1 ~> m3 s-1 or kg s-1]
    uhGM => NULL(), & !< Isopycnal height diffusion induced zonal volume fluxes [H L2 T-1 ~> m3 s-1 or kg s-1]
    vhGM => NULL()    !< Isopycnal height diffusion induced meridional volume fluxes [H L2 T-1 ~> m3 s-1 or kg s-1]

! Each of the following fields is found at nz+1 interfaces.
  real, pointer :: diapyc_vel(:,:,:) => NULL() !< The net diapycnal velocity [H s-1 ~> m s-1 or kg m-2 s-1]

end type cont_diag_ptrs

!> Vertical viscosities, drag coefficients, and related fields.
type, public :: vertvisc_type
  real :: Prandtl_turb       !< The Prandtl number for the turbulent diffusion
                             !! that is captured in Kd_shear [nondim].
  real, pointer, dimension(:,:) :: &
    bbl_thick_u => NULL(), & !< The bottom boundary layer thickness at the u-points [Z ~> m].
    bbl_thick_v => NULL(), & !< The bottom boundary layer thickness at the v-points [Z ~> m].
    kv_bbl_u => NULL(), &    !< The bottom boundary layer viscosity at the u-points [Z2 T-1 ~> m2 s-1].
    kv_bbl_v => NULL(), &    !< The bottom boundary layer viscosity at the v-points [Z2 T-1 ~> m2 s-1].
    ustar_BBL => NULL()      !< The turbulence velocity in the bottom boundary layer at h points [Z T-1 ~> m s-1].
  real, pointer, dimension(:,:) :: TKE_BBL => NULL()
                             !< A term related to the bottom boundary layer source of turbulent kinetic
                             !! energy, currently in [Z3 T-3 ~> m3 s-3], but may at some time be changed
                             !! to [R Z3 T-3 ~> W m-2].
  real, pointer, dimension(:,:) :: &
    taux_shelf => NULL(), &  !< The zonal stresses on the ocean under shelves [R Z L T-2 ~> Pa].
    tauy_shelf => NULL()     !< The meridional stresses on the ocean under shelves [R Z L T-2 ~> Pa].
  real, pointer, dimension(:,:) :: tbl_thick_shelf_u => NULL()
                !< Thickness of the viscous top boundary layer under ice shelves at u-points [Z ~> m].
  real, pointer, dimension(:,:) :: tbl_thick_shelf_v => NULL()
                !< Thickness of the viscous top boundary layer under ice shelves at v-points [Z ~> m].
  real, pointer, dimension(:,:) :: kv_tbl_shelf_u => NULL()
                !< Viscosity in the viscous top boundary layer under ice shelves at u-points [Z2 T-1 ~> m2 s-1].
  real, pointer, dimension(:,:) :: kv_tbl_shelf_v => NULL()
                !< Viscosity in the viscous top boundary layer under ice shelves at v-points [Z2 T-1 ~> m2 s-1].
  real, pointer, dimension(:,:) :: nkml_visc_u => NULL()
                !< The number of layers in the viscous surface mixed layer at u-points [nondim].
                !! This is not an integer because there may be fractional layers, and it is stored in
                !! terms of layers, not depth, to facilitate the movement of the viscous boundary layer
                !! with the flow.
  real, pointer, dimension(:,:) :: nkml_visc_v => NULL()
                !< The number of layers in the viscous surface mixed layer at v-points [nondim].
  real, pointer, dimension(:,:) :: &
    MLD => NULL()      !< Instantaneous active mixing layer depth [Z ~> m].
  real, pointer, dimension(:,:,:) :: &
    Ray_u => NULL(), & !< The Rayleigh drag velocity to be applied to each layer at u-points [Z T-1 ~> m s-1].
    Ray_v => NULL()    !< The Rayleigh drag velocity to be applied to each layer at v-points [Z T-1 ~> m s-1].
  real, pointer, dimension(:,:,:) :: Kd_shear => NULL()
                !< The shear-driven turbulent diapycnal diffusivity at the interfaces between layers
                !! in tracer columns [Z2 T-1 ~> m2 s-1].
  real, pointer, dimension(:,:,:) :: Kv_shear => NULL()
                !< The shear-driven turbulent vertical viscosity at the interfaces between layers
                !! in tracer columns [Z2 T-1 ~> m2 s-1].
  real, pointer, dimension(:,:,:) :: Kv_shear_Bu => NULL()
                !< The shear-driven turbulent vertical viscosity at the interfaces between layers in
                !! corner columns [Z2 T-1 ~> m2 s-1].
  real, pointer, dimension(:,:,:) :: Kv_slow  => NULL()
                !< The turbulent vertical viscosity component due to "slow" processes (e.g., tidal,
                !! background, convection etc) [Z2 T-1 ~> m2 s-1].
  real, pointer, dimension(:,:,:) :: TKE_turb => NULL()
                !< The turbulent kinetic energy per unit mass at the interfaces [Z2 T-2 ~> m2 s-2].
                !! This may be at the tracer or corner points
end type vertvisc_type

!> Container for information about the summed layer transports
!! and how they will vary as the barotropic velocity is changed.
type, public :: BT_cont_type
  real, allocatable :: FA_u_EE(:,:) !< The effective open face area for zonal barotropic transport
                                    !! drawing from locations far to the east [H L ~> m2 or kg m-1].
  real, allocatable :: FA_u_E0(:,:) !< The effective open face area for zonal barotropic transport
                                    !! drawing from nearby to the east [H L ~> m2 or kg m-1].
  real, allocatable :: FA_u_W0(:,:) !< The effective open face area for zonal barotropic transport
                                    !! drawing from nearby to the west [H L ~> m2 or kg m-1].
  real, allocatable :: FA_u_WW(:,:) !< The effective open face area for zonal barotropic transport
                                    !! drawing from locations far to the west [H L ~> m2 or kg m-1].
  real, allocatable :: uBT_WW(:,:)  !< uBT_WW is the barotropic velocity [L T-1 ~> m s-1], beyond which the marginal
                                    !! open face area is FA_u_WW.  uBT_WW must be non-negative.
  real, allocatable :: uBT_EE(:,:)  !< uBT_EE is a barotropic velocity [L T-1 ~> m s-1], beyond which the marginal
                                    !! open face area is FA_u_EE. uBT_EE must be non-positive.
  real, allocatable :: FA_v_NN(:,:) !< The effective open face area for meridional barotropic transport
                                    !! drawing from locations far to the north [H L ~> m2 or kg m-1].
  real, allocatable :: FA_v_N0(:,:) !< The effective open face area for meridional barotropic transport
                                    !! drawing from nearby to the north [H L ~> m2 or kg m-1].
  real, allocatable :: FA_v_S0(:,:) !< The effective open face area for meridional barotropic transport
                                    !! drawing from nearby to the south [H L ~> m2 or kg m-1].
  real, allocatable :: FA_v_SS(:,:) !< The effective open face area for meridional barotropic transport
                                    !! drawing from locations far to the south [H L ~> m2 or kg m-1].
  real, allocatable :: vBT_SS(:,:)  !< vBT_SS is the barotropic velocity, [L T-1 ~> m s-1], beyond which the marginal
                                    !! open face area is FA_v_SS. vBT_SS must be non-negative.
  real, allocatable :: vBT_NN(:,:)  !< vBT_NN is the barotropic velocity, [L T-1 ~> m s-1], beyond which the marginal
                                    !! open face area is FA_v_NN.  vBT_NN must be non-positive.
  real, allocatable :: h_u(:,:,:)   !< An effective thickness at zonal faces [H ~> m or kg m-2].
  real, allocatable :: h_v(:,:,:)   !< An effective thickness at meridional faces [H ~> m or kg m-2].
  type(group_pass_type) :: pass_polarity_BT !< Structure for polarity group halo updates
  type(group_pass_type) :: pass_FA_uv !< Structure for face area group halo updates
end type BT_cont_type

contains

!> Allocates the fields for the surface (return) properties of
!! the ocean model. Unused fields are unallocated.
subroutine allocate_surface_state(sfc_state, G, use_temperature, do_integrals, &
                                  gas_fields_ocn, use_meltpot, use_iceshelves, omit_frazil)
  type(ocean_grid_type), intent(in)    :: G                !< ocean grid structure
  type(surface),         intent(inout) :: sfc_state        !< ocean surface state type to be allocated.
  logical,     optional, intent(in)    :: use_temperature  !< If true, allocate the space for thermodynamic variables.
  logical,     optional, intent(in)    :: do_integrals     !< If true, allocate the space for vertically
                                                           !! integrated fields.
  type(coupler_1d_bc_type), &
               optional, intent(in)    :: gas_fields_ocn  !< If present, this type describes the ocean
                                              !! ocean and surface-ice fields that will participate
                                              !! in the calculation of additional gas or other
                                              !! tracer fluxes, and can be used to spawn related
                                              !! internal variables in the ice model.
  logical,     optional, intent(in)    :: use_meltpot      !< If true, allocate the space for melt potential
  logical,     optional, intent(in)    :: use_iceshelves   !< If true, allocate the space for the stresses
                                                           !! under ice shelves.
  logical,     optional, intent(in)    :: omit_frazil      !< If present and false, do not allocate the space to
                                                           !! pass frazil fluxes to the coupler

  ! local variables
  logical :: use_temp, alloc_integ, use_melt_potential, alloc_iceshelves, alloc_frazil
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: isdB, iedB, jsdB, jedB

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isdB = G%isdB ; iedB = G%iedB; jsdB = G%jsdB ; jedB = G%jedB

  use_temp = .true. ; if (present(use_temperature)) use_temp = use_temperature
  alloc_integ = .true. ; if (present(do_integrals)) alloc_integ = do_integrals
  use_melt_potential = .false. ; if (present(use_meltpot)) use_melt_potential = use_meltpot
  alloc_iceshelves = .false. ; if (present(use_iceshelves)) alloc_iceshelves = use_iceshelves
  alloc_frazil = .true. ; if (present(omit_frazil)) alloc_frazil = .not.omit_frazil

  if (sfc_state%arrays_allocated) return

  if (use_temp) then
    allocate(sfc_state%SST(isd:ied,jsd:jed)) ; sfc_state%SST(:,:) = 0.0
    allocate(sfc_state%SSS(isd:ied,jsd:jed)) ; sfc_state%SSS(:,:) = 0.0
  else
    allocate(sfc_state%sfc_density(isd:ied,jsd:jed)) ; sfc_state%sfc_density(:,:) = 0.0
  endif
  if (use_temp .and. alloc_frazil) then
    allocate(sfc_state%frazil(isd:ied,jsd:jed)) ; sfc_state%frazil(:,:) = 0.0
  endif
  allocate(sfc_state%sea_lev(isd:ied,jsd:jed)) ; sfc_state%sea_lev(:,:) = 0.0
  allocate(sfc_state%Hml(isd:ied,jsd:jed)) ; sfc_state%Hml(:,:) = 0.0
  allocate(sfc_state%u(IsdB:IedB,jsd:jed)) ; sfc_state%u(:,:) = 0.0
  allocate(sfc_state%v(isd:ied,JsdB:JedB)) ; sfc_state%v(:,:) = 0.0

  if (use_melt_potential) then
    allocate(sfc_state%melt_potential(isd:ied,jsd:jed)) ; sfc_state%melt_potential(:,:) = 0.0
  endif

  if (alloc_integ) then
    ! Allocate structures for the vertically integrated ocean_mass, ocean_heat, and ocean_salt.
    allocate(sfc_state%ocean_mass(isd:ied,jsd:jed)) ; sfc_state%ocean_mass(:,:) = 0.0
    if (use_temp) then
      allocate(sfc_state%ocean_heat(isd:ied,jsd:jed)) ; sfc_state%ocean_heat(:,:) = 0.0
      allocate(sfc_state%ocean_salt(isd:ied,jsd:jed)) ; sfc_state%ocean_salt(:,:) = 0.0
      allocate(sfc_state%TempxPmE(isd:ied,jsd:jed))   ; sfc_state%TempxPmE(:,:) = 0.0
      allocate(sfc_state%salt_deficit(isd:ied,jsd:jed))  ; sfc_state%salt_deficit(:,:) = 0.0
      allocate(sfc_state%internal_heat(isd:ied,jsd:jed)) ; sfc_state%internal_heat(:,:) = 0.0
    endif
  endif

  if (alloc_iceshelves) then
    allocate(sfc_state%taux_shelf(IsdB:IedB,jsd:jed)) ; sfc_state%taux_shelf(:,:) = 0.0
    allocate(sfc_state%tauy_shelf(isd:ied,JsdB:JedB)) ; sfc_state%tauy_shelf(:,:) = 0.0
  endif

  if (present(gas_fields_ocn)) &
    call coupler_type_spawn(gas_fields_ocn, sfc_state%tr_fields, &
                            (/is,is,ie,ie/), (/js,js,je,je/), as_needed=.true.)

  sfc_state%arrays_allocated = .true.

end subroutine allocate_surface_state

!> Deallocates the elements of a surface state type.
subroutine deallocate_surface_state(sfc_state)
  type(surface), intent(inout) :: sfc_state !< ocean surface state type to be deallocated here.

  if (.not.sfc_state%arrays_allocated) return

  if (allocated(sfc_state%melt_potential)) deallocate(sfc_state%melt_potential)
  if (allocated(sfc_state%SST)) deallocate(sfc_state%SST)
  if (allocated(sfc_state%SSS)) deallocate(sfc_state%SSS)
  if (allocated(sfc_state%sfc_density)) deallocate(sfc_state%sfc_density)
  if (allocated(sfc_state%sea_lev)) deallocate(sfc_state%sea_lev)
  if (allocated(sfc_state%Hml)) deallocate(sfc_state%Hml)
  if (allocated(sfc_state%u)) deallocate(sfc_state%u)
  if (allocated(sfc_state%v)) deallocate(sfc_state%v)
  if (allocated(sfc_state%ocean_mass)) deallocate(sfc_state%ocean_mass)
  if (allocated(sfc_state%ocean_heat)) deallocate(sfc_state%ocean_heat)
  if (allocated(sfc_state%ocean_salt)) deallocate(sfc_state%ocean_salt)
  if (allocated(sfc_state%salt_deficit)) deallocate(sfc_state%salt_deficit)

  call coupler_type_destructor(sfc_state%tr_fields)

  sfc_state%arrays_allocated = .false.

end subroutine deallocate_surface_state

!> Rotate the surface state fields from the input to the model indices.
subroutine rotate_surface_state(sfc_state_in, G_in, sfc_state, G, turns)
  type(surface), intent(in) :: sfc_state_in
  type(ocean_grid_type), intent(in) :: G_in
  type(surface), intent(inout) :: sfc_state
  type(ocean_grid_type), intent(in) :: G
  integer, intent(in) :: turns

  logical :: use_temperature, do_integrals, use_melt_potential, use_iceshelves

  ! NOTE: Many of these are weak tests, since only one is checked
  use_temperature = allocated(sfc_state_in%SST) &
      .and. allocated(sfc_state_in%SSS)
  use_melt_potential = allocated(sfc_state_in%melt_potential)
  do_integrals = allocated(sfc_state_in%ocean_mass)
  use_iceshelves = allocated(sfc_state_in%taux_shelf) &
      .and. allocated(sfc_state_in%tauy_shelf)

  if (.not. sfc_state%arrays_allocated) then
    call allocate_surface_state(sfc_state, G, &
        use_temperature=use_temperature, &
        do_integrals=do_integrals, &
        use_meltpot=use_melt_potential, &
        use_iceshelves=use_iceshelves &
    )
    sfc_state%arrays_allocated = .true.
  endif

  if (use_temperature) then
    call rotate_array(sfc_state_in%SST, turns, sfc_state%SST)
    call rotate_array(sfc_state_in%SSS, turns, sfc_state%SSS)
  else
    call rotate_array(sfc_state_in%sfc_density, turns, sfc_state%sfc_density)
  endif

  call rotate_array(sfc_state_in%Hml, turns, sfc_state%Hml)
  call rotate_vector(sfc_state_in%u, sfc_state_in%v, turns, &
      sfc_state%u, sfc_state%v)
  call rotate_array(sfc_state_in%sea_lev, turns, sfc_state%sea_lev)

  if (use_melt_potential) then
    call rotate_array(sfc_state_in%melt_potential, turns, sfc_state%melt_potential)
  endif

  if (do_integrals) then
    call rotate_array(sfc_state_in%ocean_mass, turns, sfc_state%ocean_mass)
    if (use_temperature) then
      call rotate_array(sfc_state_in%ocean_heat, turns, sfc_state%ocean_heat)
      call rotate_array(sfc_state_in%ocean_salt, turns, sfc_state%ocean_salt)
      call rotate_array(sfc_state_in%SSS, turns, sfc_state%SSS)
      call rotate_array(sfc_state_in%salt_deficit, turns, sfc_state%salt_deficit)
      call rotate_array(sfc_state_in%internal_heat, turns, sfc_state%internal_heat)
    endif
  endif

  if (use_iceshelves) then
    call rotate_vector(sfc_state_in%taux_shelf, sfc_state_in%tauy_shelf, turns, &
        sfc_state%taux_shelf, sfc_state%tauy_shelf)
  endif

  if (use_temperature .and. allocated(sfc_state_in%frazil)) &
    call rotate_array(sfc_state_in%frazil, turns, sfc_state%frazil)

  ! Scalar transfers
  sfc_state%T_is_conT = sfc_state_in%T_is_conT
  sfc_state%S_is_absS = sfc_state_in%S_is_absS

  ! TODO: tracer field rotation
  if (coupler_type_initialized(sfc_state_in%tr_fields)) &
    call MOM_error(FATAL, "Rotation of surface state tracers is not yet implemented.")
end subroutine rotate_surface_state

!> Allocates the arrays contained within a BT_cont_type and initializes them to 0.
subroutine alloc_BT_cont_type(BT_cont, G, GV, alloc_faces)
  type(BT_cont_type),      pointer    :: BT_cont !< The BT_cont_type whose elements will be allocated
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure.
  logical,       optional, intent(in) :: alloc_faces !< If present and true, allocate
                                            !! memory for effective face thicknesses.

  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(BT_cont)) call MOM_error(FATAL, &
    "alloc_BT_cont_type called with an associated BT_cont_type pointer.")

  allocate(BT_cont)
  allocate(BT_cont%FA_u_WW(IsdB:IedB,jsd:jed)) ; BT_cont%FA_u_WW(:,:) = 0.0
  allocate(BT_cont%FA_u_W0(IsdB:IedB,jsd:jed)) ; BT_cont%FA_u_W0(:,:) = 0.0
  allocate(BT_cont%FA_u_E0(IsdB:IedB,jsd:jed)) ; BT_cont%FA_u_E0(:,:) = 0.0
  allocate(BT_cont%FA_u_EE(IsdB:IedB,jsd:jed)) ; BT_cont%FA_u_EE(:,:) = 0.0
  allocate(BT_cont%uBT_WW(IsdB:IedB,jsd:jed))  ; BT_cont%uBT_WW(:,:) = 0.0
  allocate(BT_cont%uBT_EE(IsdB:IedB,jsd:jed))  ; BT_cont%uBT_EE(:,:) = 0.0

  allocate(BT_cont%FA_v_SS(isd:ied,JsdB:JedB)) ; BT_cont%FA_v_SS(:,:) = 0.0
  allocate(BT_cont%FA_v_S0(isd:ied,JsdB:JedB)) ; BT_cont%FA_v_S0(:,:) = 0.0
  allocate(BT_cont%FA_v_N0(isd:ied,JsdB:JedB)) ; BT_cont%FA_v_N0(:,:) = 0.0
  allocate(BT_cont%FA_v_NN(isd:ied,JsdB:JedB)) ; BT_cont%FA_v_NN(:,:) = 0.0
  allocate(BT_cont%vBT_SS(isd:ied,JsdB:JedB))  ; BT_cont%vBT_SS(:,:) = 0.0
  allocate(BT_cont%vBT_NN(isd:ied,JsdB:JedB))  ; BT_cont%vBT_NN(:,:) = 0.0

  if (present(alloc_faces)) then ; if (alloc_faces) then
    allocate(BT_cont%h_u(IsdB:IedB,jsd:jed,1:nz)) ; BT_cont%h_u(:,:,:) = 0.0
    allocate(BT_cont%h_v(isd:ied,JsdB:JedB,1:nz)) ; BT_cont%h_v(:,:,:) = 0.0
  endif ; endif

end subroutine alloc_BT_cont_type

!> Deallocates the arrays contained within a BT_cont_type.
subroutine dealloc_BT_cont_type(BT_cont)
  type(BT_cont_type), pointer :: BT_cont !< The BT_cont_type whose elements will be deallocated.

  if (.not.associated(BT_cont)) return

  deallocate(BT_cont%FA_u_WW) ; deallocate(BT_cont%FA_u_W0)
  deallocate(BT_cont%FA_u_E0) ; deallocate(BT_cont%FA_u_EE)
  deallocate(BT_cont%uBT_WW)  ; deallocate(BT_cont%uBT_EE)

  deallocate(BT_cont%FA_v_SS) ; deallocate(BT_cont%FA_v_S0)
  deallocate(BT_cont%FA_v_N0) ; deallocate(BT_cont%FA_v_NN)
  deallocate(BT_cont%vBT_SS)  ; deallocate(BT_cont%vBT_NN)

  if (allocated(BT_cont%h_u)) deallocate(BT_cont%h_u)
  if (allocated(BT_cont%h_v)) deallocate(BT_cont%h_v)

  deallocate(BT_cont)

end subroutine dealloc_BT_cont_type

!> Diagnostic checksums on various elements of a thermo_var_ptrs type for debugging.
subroutine MOM_thermovar_chksum(mesg, tv, G)
  character(len=*),      intent(in) :: mesg !< A message that appears in the checksum lines
  type(thermo_var_ptrs), intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  if (associated(tv%T)) &
    call hchksum(tv%T, mesg//" tv%T", G%HI)
  if (associated(tv%S)) &
    call hchksum(tv%S, mesg//" tv%S", G%HI)
  if (associated(tv%frazil)) &
    call hchksum(tv%frazil, mesg//" tv%frazil", G%HI, scale=G%US%Q_to_J_kg*G%US%RZ_to_kg_m2)
  if (associated(tv%salt_deficit)) &
    call hchksum(tv%salt_deficit, mesg//" tv%salt_deficit", G%HI, scale=G%US%RZ_to_kg_m2)
  if (associated(tv%TempxPmE)) &
    call hchksum(tv%TempxPmE, mesg//" tv%TempxPmE", G%HI, scale=G%US%RZ_to_kg_m2)
end subroutine MOM_thermovar_chksum

end module MOM_variables
