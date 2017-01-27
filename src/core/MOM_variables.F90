module MOM_variables

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

use MOM_domains, only : MOM_domain_type, get_domain_extent, group_pass_type
use MOM_debugging, only : hchksum
use MOM_error_handler, only : MOM_error, FATAL
use MOM_grid, only : ocean_grid_type
use MOM_io, only : vardesc
use MOM_EOS, only : EOS_type

use coupler_types_mod, only : coupler_2d_bc_type

implicit none ; private

#include <MOM_memory.h>

public MOM_thermovar_chksum
public ocean_grid_type, vardesc, alloc_BT_cont_type, dealloc_BT_cont_type

type, public :: p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d
type, public :: p2d
  real, dimension(:,:), pointer :: p => NULL()
end type p2d

!   The following structure contains pointers to various fields
! which may be used describe the surface state of MOM, and which
! will be returned to a the calling program
type, public :: surface
  real, pointer, dimension(:,:) :: &
    SST => NULL(), &     ! The sea surface temperature in C.
    SSS => NULL(), &     ! The sea surface salinity in psu.
    sfc_density => NULL(), &  ! The mixed layer density in kg m-3.
    u => NULL(), &       ! The mixed layer zonal velocity in m s-1.
    v => NULL(), &       ! The mixed layer meridional velocity in m s-1.
    Hml => NULL(), &     ! The mixed layer depth in m.
    ocean_mass => NULL(), &  ! The total mass of the ocean in kg m-2.
    ocean_heat => NULL(), &  ! The total heat content of the ocean in C kg m-2.
    ocean_salt => NULL(), &  ! The total salt content of the ocean in kgSalt m-2.
    taux_shelf => NULL(), &  ! The zonal and meridional stresses on the ocean
    tauy_shelf => NULL(), &  ! under shelves, in Pa.
    frazil => NULL(), &  ! The energy needed to heat the ocean column to the
                         ! freezing point over the call to step_MOM, in J m-2.
    salt_deficit => NULL(), & ! The salt needed to maintain the ocean column
                         ! at a minimum salinity of 0.01 PSU over the call to
                         ! step_MOM, in kgSalt m-2.
    TempxPmE => NULL(), &  ! The net inflow of water into the ocean times
                         ! the temperature at which this inflow occurs during
                         ! the call to step_MOM, in deg C kg m-2.
                         !   This should be prescribed in the forcing fields,
                         ! but as it often is not, this is a useful heat budget
                         ! diagnostic.
    internal_heat => NULL() , & ! Any internal or geothermal heat sources that
                         ! are applied to the ocean integrated over the call
                         ! to step_MOM, in deg C kg m-2.
    sea_lev => NULL()    ! The sea level in m.  If a reduced  surface gravity is
                         ! used, that is compensated for in sea_lev.
  type(coupler_2d_bc_type), pointer :: tr_fields  => NULL()
                                          ! A structure that may contain an
                                          ! array of named fields describing
                                          ! tracer-related quantities.
       !!! NOTE: ALL OF THE ARRAYS IN TR_FIELDS USE THE COUPLER'S INDEXING
       !!!       CONVENTION AND HAVE NO HALOS!  THIS IS DONE TO CONFORM TO
       !!!       THE TREATMENT IN MOM4, BUT I DON'T LIKE IT!
end type surface

!   The following structure contains pointers to an assortment of
! thermodynamic fields that may be available, including potential
! temperature, salinity and mixed layer density.
type, public :: thermo_var_ptrs
!   If allocated, the following variables have nz layers.
  real, pointer :: T(:,:,:) => NULL()   ! Potential temperature in C.
  real, pointer :: S(:,:,:) => NULL()   ! Salnity in psu.
  type(EOS_type), pointer :: eqn_of_state => NULL() ! Type that indicates the
                                        ! equation of state to use.
!  integer :: nk_Rml = 0  !   The number of variable density near-surface layers
!                         ! with a bulk mixed layer (nkml+nkbl).
!  integer :: nkml = 0    !   The number of layers within a surface bulk mixed
!                         ! layer.
  real :: P_Ref          !   The coordinate-density reference pressure in Pa.
                         ! This is the pressure used to calculate Rml from
                         ! T and S when eqn_of_state is associated.
  real :: C_p            !   The heat capacity of seawater, in J K-1 kg-1.
                         ! When conservative temperature is used, this is
                         ! constant and exactly 3991.86795711963 J K kg-1.
  real, pointer, dimension(:,:) :: &
    Hml => NULL(), &     !   The surface mixed layer depth in m.

!  These arrays are accumulated fluxes for communication with other components.
    frazil => NULL(), &  !   The energy needed to heat the ocean column to the
                         ! freezing point since calculate_surface_state was
                         ! last called, in units of J m-2.
    salt_deficit => NULL(), & !   The salt needed to maintain the ocean column
                         ! at a minumum salinity of 0.01 PSU since the last time
                         ! that calculate_surface_state was called, in units
                         ! of gSalt m-2.
    TempxPmE => NULL(), &!   The net inflow of water into the ocean times the
                         ! temperature at which this inflow occurs since the
                         ! last call to calculate_surface_state, in units of
                         ! deg C kg m-2. This should be prescribed in the
                         ! forcing fields, but as it often is not, this is a
                         ! useful heat budget diagnostic.
    internal_heat => NULL() ! Any internal or geothermal heat sources that
                         ! have been applied to the ocean since the last call to
                         ! calculate_surface_state, in units of deg C kg m-2.
end type thermo_var_ptrs

type, public :: ocean_internal_state
! This structure contains pointers to all of the prognostic variables allocated
! here and in MOM.F90.  It is useful for sending them for diagnostics, and for
! preparation for ensembles later on.  All variables have the same names as the
! local (public) variables they refer to in MOM.F90.
  real, pointer, dimension(:,:,:) :: &
    u => NULL(), v => NULL(), h => NULL()
  real, pointer, dimension(:,:,:) :: &
    uh => NULL(), vh => NULL(), &
    CAu => NULL(), CAv => NULL(), &
    PFu  => NULL(), PFv => NULL(), diffu => NULL(), diffv => NULL(), &
    T => NULL(), S => NULL(), &
    pbce => NULL(), u_accel_bt => NULL(), v_accel_bt => NULL(), &
    u_av => NULL(), v_av => NULL(), u_prev => NULL(), v_prev => NULL()
end type ocean_internal_state

type, public :: accel_diag_ptrs
! This structure contains pointers to arrays with accelerations, which can
! later be used for derived diagnostics, like energy balances.

! Each of the following fields has nz layers.
  real, pointer :: diffu(:,:,:) => NULL()    ! Accelerations due to along iso-
  real, pointer :: diffv(:,:,:) => NULL()    ! pycnal viscosity, in m s-2.
  real, pointer :: CAu(:,:,:) => NULL()      ! Coriolis and momentum advection
  real, pointer :: CAv(:,:,:) => NULL()      ! accelerations, in m s-2.
  real, pointer :: PFu(:,:,:) => NULL()      ! Accelerations due to pressure
  real, pointer :: PFv(:,:,:) => NULL()      ! forces, in m s-2.
  real, pointer :: du_dt_visc(:,:,:) => NULL()! Accelerations due to vertical
  real, pointer :: dv_dt_visc(:,:,:) => NULL()! viscosity, in m s-2.
  real, pointer :: du_dt_dia(:,:,:) => NULL()! Accelerations due to diapycnal
  real, pointer :: dv_dt_dia(:,:,:) => NULL()! mixing, in m s-2.
  real, pointer :: du_other(:,:,:) => NULL() ! Velocity changes due to any other
  real, pointer :: dv_other(:,:,:) => NULL() ! processes that are not due to any
                                             ! explicit accelerations, in m s-1.

  ! These accelerations are sub-terms included in the accelerations above.
  real, pointer :: gradKEu(:,:,:) => NULL()  ! gradKEu = - d/dx(u2), in m s-2.
  real, pointer :: gradKEv(:,:,:) => NULL()  ! gradKEv = - d/dy(u2), in m s-2.
  real, pointer :: rv_x_v(:,:,:) => NULL()   ! rv_x_v = rv * v at u, in m s-2.
  real, pointer :: rv_x_u(:,:,:) => NULL()   ! rv_x_u = rv * u at v, in m s-2.

end type accel_diag_ptrs

type, public :: cont_diag_ptrs
! This structure contains pointers to arrays with accelerations, which can
! later be used for derived diagnostics, like energy balances.

! Each of the following fields has nz layers.
  real, pointer :: uh(:,:,:) => NULL()    ! Resolved layer thickness fluxes,
  real, pointer :: vh(:,:,:) => NULL()    ! in m3 s-1 or kg s-1.
  real, pointer :: uhGM(:,:,:) => NULL()  ! Thickness diffusion induced
  real, pointer :: vhGM(:,:,:) => NULL()  ! volume fluxes in m3 s-1.

! Each of the following fields is found at nz+1 interfaces.
  real, pointer :: diapyc_vel(:,:,:) => NULL()! The net diapycnal velocity,

end type cont_diag_ptrs

type, public :: vertvisc_type
!   This structure contains vertical viscosities, drag coefficients, and
! related fields.
  logical :: calc_bbl            ! If true, the BBL viscosity and thickness
                                 ! need to be recalculated.
  real :: bbl_calc_time_interval ! The amount of time over which the impending
                                 ! calculation of the BBL properties will apply,
                                 ! for use in diagnostics of the BBL properties.
  real :: Prandtl_turb       ! The Prandtl number for the turbulent diffusion
                             ! that is captured in Kd_turb.
  real, pointer, dimension(:,:) :: &
    bbl_thick_u => NULL(), & ! The bottom boundary layer thickness at the zonal
    bbl_thick_v => NULL(), & ! and meridional velocity points, in m.
    kv_bbl_u => NULL(), &    ! The bottom boundary layer viscosity at the zonal
    kv_bbl_v => NULL(), &    ! and meridional velocity points, in m2 s-1.
    ustar_BBL => NULL(), &   ! The turbulence velocity in the bottom boundary
                             ! layer at h points, in m s-1.
    TKE_BBL => NULL(), &     ! A term related to the bottom boundary layer
                             ! source of turbulent kinetic energy, currently
                             ! in units of m3 s-3, but will later be changed
                             ! to W m-2.
    taux_shelf => NULL(), &  ! The zonal and meridional stresses on the ocean
    tauy_shelf => NULL(), &  ! under shelves, in Pa.
    tbl_thick_shelf_u => NULL(), & ! Thickness of the viscous top boundary
    tbl_thick_shelf_v => NULL(), & ! layer under ice shelves, in m.
    kv_tbl_shelf_u => NULL(), &  ! Viscosity in the viscous top boundary
    kv_tbl_shelf_v => NULL(), &  ! layer under ice shelves, in m2 s-1.
    nkml_visc_u => NULL(), & ! The number of layers in the viscous surface
    nkml_visc_v => NULL(), & ! mixed layer.  These are not integers because
                             ! there may be fractional layers.  This is done
                             ! in terms of layers, not depth, to facilitate
                             ! the movement of the viscous boundary layer with
                             ! the flow.
    MLD => NULL()            !< Instantaneous active mixing layer depth (H units).
  real, pointer, dimension(:,:,:) :: &
    Ray_u => NULL(), &  ! The Rayleigh drag velocity to be applied to each layer
    Ray_v => NULL(), &  ! at u- and v-points, in m s-1.
    Kd_extra_T => NULL(), & ! The extra diffusivities of temperature and
    Kd_extra_S => NULL(), & ! salinity due to double diffusion at the interfaces
                        ! relative to the diffusivity of density, in m2 s-1.
                        ! One of these is always 0.  Kd_extra_S is positive for
                        ! salt fingering; Kd_extra_T is positive for double
                        ! diffusive convection.  These are only allocated if
                        ! DOUBLE_DIFFUSION is true.
    Kd_turb => NULL(), &! The turbulent diapycnal diffusivity at the interfaces
                        ! between each layer, in m2 s-1.
    Kv_turb => NULL(), &! The turbulent vertical viscosity at the interfaces
                        ! between each layer, in m2 s-1.
    TKE_turb => NULL()  ! The turbulent kinetic energy per unit mass defined
                        ! at the interfaces between each layer, in m2 s-2.
end type vertvisc_type

type, public :: BT_cont_type
  real, pointer, dimension(:,:) :: &
    FA_u_EE => NULL(), &  ! The FA_u_XX variables are the effective open face
    FA_u_E0 => NULL(), &  ! areas for barotropic transport through the zonal
    FA_u_W0 => NULL(), &  ! faces, all in H m, with the XX indicating where
    FA_u_WW => NULL(), &  ! the transport is from, with _EE drawing from points
                          ! far to the east, _E0 from points nearby from the
                          ! east, _W0 nearby from the west, and _WW from far to
                          ! the west.
    uBT_WW => NULL(), &   ! uBT_WW is the barotropic velocity, in m s-1, beyond
                          ! which the marginal open face area is FA_u_WW.
                          ! uBT_EE must be non-negative.
    uBT_EE => NULL(), &   ! uBT_EE is the barotropic velocity, in m s-1, beyond
                          ! which the marginal open face area is FA_u_EE.
                          ! uBT_EE must be non-positive.
    FA_v_NN => NULL(), &  ! The FA_v_XX variables are the effective open face
    FA_v_N0 => NULL(), &  ! areas for barotropic transport through the meridional
    FA_v_S0 => NULL(), &  ! faces, all in H m, with the XX indicating where
    FA_v_SS => NULL(), &  ! the transport is from, with _NN drawing from points
                          ! far to the north, _N0 from points nearby from the
                          ! north, _S0 nearby from the south, and _SS from far
                          ! to the south.
    vBT_SS => NULL(), &   ! vBT_SS is the barotropic velocity, in m s-1, beyond
                          ! which the marginal open face area is FA_v_SS.
                          ! vBT_SS must be non-negative.
    vBT_NN => NULL()      ! vBT_NN is the barotropic velocity, in m s-1, beyond
                          ! which the marginal open face area is FA_v_NN.
                          ! vBT_NN must be non-positive.
  real, pointer, dimension(:,:,:) :: &
    h_u => NULL(), &      ! An effective thickness at zonal faces, in H.
    h_v => NULL()         ! An effective thickness at meridional faces, in H.
  type(group_pass_type) :: pass_polarity_BT, pass_FA_uv ! For group halo updates
end type BT_cont_type

contains


subroutine alloc_BT_cont_type(BT_cont, G, alloc_faces)
  type(BT_cont_type),    pointer    :: BT_cont
  type(ocean_grid_type), intent(in) :: G
  logical,     optional, intent(in) :: alloc_faces

  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
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
    allocate(BT_cont%h_u(IsdB:IedB,jsd:jed,1:G%ke)) ; BT_cont%h_u(:,:,:) = 0.0
    allocate(BT_cont%h_v(isd:ied,JsdB:JedB,1:G%ke)) ; BT_cont%h_v(:,:,:) = 0.0
  endif ; endif

end subroutine alloc_BT_cont_type

subroutine dealloc_BT_cont_type(BT_cont)
  type(BT_cont_type), pointer :: BT_cont

  if (.not.associated(BT_cont)) return

  deallocate(BT_cont%FA_u_WW) ; deallocate(BT_cont%FA_u_W0)
  deallocate(BT_cont%FA_u_E0) ; deallocate(BT_cont%FA_u_EE)
  deallocate(BT_cont%uBT_WW)  ; deallocate(BT_cont%uBT_EE)

  deallocate(BT_cont%FA_v_SS) ; deallocate(BT_cont%FA_v_S0)
  deallocate(BT_cont%FA_v_N0) ; deallocate(BT_cont%FA_v_NN)
  deallocate(BT_cont%vBT_SS)  ; deallocate(BT_cont%vBT_NN)

  if (associated(BT_cont%h_u)) deallocate(BT_cont%h_u)
  if (associated(BT_cont%h_v)) deallocate(BT_cont%h_v)

  deallocate(BT_cont)

end subroutine dealloc_BT_cont_type

subroutine MOM_thermovar_chksum(mesg, tv, G)
  character(len=*),                    intent(in) :: mesg
  type(thermo_var_ptrs),               intent(in) :: tv
  type(ocean_grid_type),               intent(in) :: G
!   This subroutine writes out chksums for the model's basic state variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m.
!  (in)      uh - Volume flux through zonal faces = u*h*dy, m3 s-1.
!  (in)      vh - Volume flux through meridional faces = v*h*dx, in m3 s-1.
!  (in)      G - The ocean's grid structure.
  integer :: is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  if (associated(tv%T)) &
    call hchksum(tv%T, mesg//" tv%T",G%HI)
  if (associated(tv%S)) &
    call hchksum(tv%S, mesg//" tv%S",G%HI)
  if (associated(tv%frazil)) &
    call hchksum(tv%frazil, mesg//" tv%frazil",G%HI)
  if (associated(tv%salt_deficit)) &
    call hchksum(tv%salt_deficit, mesg//" tv%salt_deficit",G%HI)
  if (associated(tv%TempxPmE)) &
    call hchksum(tv%TempxPmE, mesg//" tv%TempxPmE",G%HI)
end subroutine MOM_thermovar_chksum

end module MOM_variables
