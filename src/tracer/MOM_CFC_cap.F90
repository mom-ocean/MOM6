!> Simulates CFCs using atmospheric pressure, wind speed and sea ice cover
!! provided via cap (only NUOPC cap is implemented so far).
module MOM_CFC_cap

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms,            only : EFP_type
use MOM_diag_mediator,   only : diag_ctrl, register_diag_field, post_data
use MOM_error_handler,   only : MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_hor_index,       only : hor_index_type
use MOM_grid,            only : ocean_grid_type
use MOM_io,              only : file_exists, MOM_read_data, slasher
use MOM_io,              only : vardesc, var_desc, query_vardesc, stdout
use MOM_open_boundary,   only : ocean_OBC_type
use MOM_restart,         only : query_initialized, MOM_restart_CS
use MOM_spatial_means,   only : global_mass_int_EFP
use MOM_time_manager,    only : time_type
use time_interp_external_mod, only : init_external_field, time_interp_external
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init,   only : tracer_Z_init
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_CFC_cap, initialize_CFC_cap, CFC_cap_unit_tests
public CFC_cap_column_physics, CFC_cap_surface_state, CFC_cap_fluxes
public CFC_cap_stock, CFC_cap_end

integer, parameter :: NTR = 2 !< the number of tracers in this module.

!> The control structure for the CFC_cap tracer package
type, public :: CFC_cap_CS ; private
  character(len=200) :: IC_file !< The file in which the CFC initial values can
                                !! be found, or an empty string for internal initilaization.
  logical :: Z_IC_file !< If true, the IC_file is in Z-space.  The default is false.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the MOM6 tracer registry
  real, pointer, dimension(:,:,:) :: &
    CFC11 => NULL(), &     !< The CFC11 concentration [mol kg-1].
    CFC12 => NULL()        !< The CFC12 concentration [mol kg-1].
  ! In the following variables a suffix of _11 refers to CFC11 and _12 to CFC12.
  real :: CFC11_IC_val = 0.0    !< The initial value assigned to CFC11 [mol kg-1].
  real :: CFC12_IC_val = 0.0    !< The initial value assigned to CFC12 [mol kg-1].
  real :: CFC11_land_val = -1.0 !< The value of CFC11 used where land is masked out [mol kg-1].
  real :: CFC12_land_val = -1.0 !< The value of CFC12 used where land is masked out [mol kg-1].
  logical :: tracers_may_reinit !< If true, tracers may be reset via the initialization code
                                !! if they are not found in the restart files.
  character(len=16) :: CFC11_name !< CFC11 variable name
  character(len=16) :: CFC12_name !< CFC12 variable name
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate
                                             !! the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< Model restart control structure

  ! The following vardesc types contain a package of metadata about each tracer.
  type(vardesc) :: CFC11_desc !< A set of metadata for the CFC11 tracer
  type(vardesc) :: CFC12_desc !< A set of metadata for the CFC12 tracer
  !>@{ Diagnostic IDs
  integer :: id_cfc11_cmor = -1, id_cfc12_cmor = -1
  !>@}
end type CFC_cap_CS

contains

!> Register the CFCs to be used with MOM and read the parameters
!! that are used with this tracer package
function register_CFC_cap(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),    intent(in) :: HI         !< A horizontal index type structure.
  type(verticalGrid_type), intent(in) :: GV         !< The ocean's vertical grid structure.
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters.
  type(CFC_cap_CS),        pointer    :: CS         !< A pointer that is set to point to the control
                                                    !! structure for this module.
  type(tracer_registry_type), &
                           pointer    :: tr_Reg     !< A pointer to the tracer registry.
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control struct

  ! Local variables
  character(len=40)  :: mdl = "MOM_CFC_cap" ! This module's name.
  character(len=200) :: inputdir ! The directory where NetCDF input files are.
  ! This include declares and sets the variable "version".
#include "version_variable.h"
  real, dimension(:,:,:), pointer :: tr_ptr => NULL()
  character(len=200) :: dummy      ! Dummy variable to store params that need to be logged here.
  logical :: register_CFC_cap
  integer :: isd, ied, jsd, jed, nz, m

  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "register_CFC_cap called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "CFC_IC_FILE", CS%IC_file, &
                 "The file in which the CFC initial values can be "//&
                 "found, or an empty string for internal initialization.", &
                 default=" ")
  if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file,'/') == 0)) then
    ! Add the directory if CS%IC_file is not already a complete path.
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    call log_param(param_file, mdl, "INPUTDIR/CFC_IC_FILE", CS%IC_file)
  endif
  call get_param(param_file, mdl, "CFC_IC_FILE_IS_Z", CS%Z_IC_file, &
                 "If true, CFC_IC_FILE is in depth space, not layer space", &
                 default=.false.)
  call get_param(param_file, mdl, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code "//&
                 "if they are not found in the restart files.  Otherwise "//&
                 "it is a fatal error if tracers are not found in the "//&
                 "restart files of a restarted run.", default=.false.)
  call get_param(param_file, mdl, "CFC11_IC_VAL", CS%CFC11_IC_val, &
                 "Value that CFC_11 is set to when it is not read from a file.", &
                 units="mol kg-1", default=0.0)
  call get_param(param_file, mdl, "CFC12_IC_VAL", CS%CFC12_IC_val, &
                 "Value that CFC_12 is set to when it is not read from a file.", &
                 units="mol kg-1", default=0.0)

  ! the following params are not used in this module. Instead, they are used in
  ! the cap but are logged here to keep all the CFC cap params together.
  call get_param(param_file, mdl, "CFC_BC_FILE", dummy, &
                "The file in which the CFC-11 and CFC-12 atm concentrations can be "//&
                "found (units must be parts per trillion), or an empty string for "//&
                "internal BC generation (TODO).", default=" ")
  if ((len_trim(dummy) > 0) .and. (scan(dummy,'/') == 0)) then
    ! Add the directory if dummy is not already a complete path.
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    dummy = trim(slasher(inputdir))//trim(dummy)
    call log_param(param_file, mdl, "INPUTDIR/CFC_IC_FILE", dummy)
  endif
  if (len_trim(dummy) > 0) then
    call get_param(param_file, mdl, "CFC11_VARIABLE", dummy, &
                 "The name of the variable representing CFC-11 in  "//&
                 "CFC_BC_FILE.", default="CFC_11")
    call get_param(param_file, mdl, "CFC12_VARIABLE", dummy, &
                 "The name of the variable representing CFC-12 in  "//&
                 "CFC_BC_FILE.", default="CFC_12")
  endif

  ! The following vardesc types contain a package of metadata about each tracer,
  ! including, the name; units; longname; and grid information.
  CS%CFC11_name = "CFC_11" ; CS%CFC12_name = "CFC_12"
  CS%CFC11_desc = var_desc(CS%CFC11_name,"mol kg-1","Moles Per Unit Mass of CFC-11 in sea water", caller=mdl)
  CS%CFC12_desc = var_desc(CS%CFC12_name,"mol kg-1","Moles Per Unit Mass of CFC-12 in sea water", caller=mdl)

  allocate(CS%CFC11(isd:ied,jsd:jed,nz), source=0.0)
  allocate(CS%CFC12(isd:ied,jsd:jed,nz), source=0.0)

  ! This pointer assignment is needed to force the compiler not to do a copy in
  ! the registration calls.  Curses on the designers and implementers of F90.
  tr_ptr => CS%CFC11
  ! Register CFC11 for horizontal advection, diffusion, and restarts.
  call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                       tr_desc=CS%CFC11_desc, registry_diags=.true., &
                       restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit)
  ! Do the same for CFC12
  tr_ptr => CS%CFC12
  call register_tracer(tr_ptr, Tr_Reg, param_file, HI, GV, &
                       tr_desc=CS%CFC12_desc, registry_diags=.true., &
                       restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit)

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_CFC_cap = .true.

end function register_CFC_cap

!> Initialize the CFC tracer fields and set up the tracer output.
subroutine initialize_CFC_cap(restart, day, G, GV, US, h, diag, OBC, CS)
  logical,                        intent(in) :: restart    !< .true. if the fields have already been
                                                           !! read from a restart file.
  type(time_type), target,        intent(in) :: day        !< Time of the start of the run.
  type(ocean_grid_type),          intent(in) :: G          !< The ocean's grid structure.
  type(verticalGrid_type),        intent(in) :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),          intent(in) :: US         !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                  intent(in) :: h          !< Layer thicknesses [H ~> m or kg m-2].
  type(diag_ctrl), target,        intent(in) :: diag       !< A structure that is used to regulate
                                                           !! diagnostic output.
  type(ocean_OBC_type),           pointer    :: OBC        !< This open boundary condition type
                                                           !! specifies whether, where, and what
                                                           !! open boundary conditions are used.
  type(CFC_cap_CS),              pointer    :: CS          !< The control structure returned by a
                                                           !! previous call to register_CFC_cap.

  ! local variables
  logical :: from_file = .false.

  if (.not.associated(CS)) return

  CS%Time => day
  CS%diag => diag

  if (.not.restart .or. (CS%tracers_may_reinit .and. &
      .not.query_initialized(CS%CFC11, CS%CFC11_name, CS%restart_CSp))) &
    call init_tracer_CFC(h, CS%CFC11, CS%CFC11_name, CS%CFC11_land_val, &
                         CS%CFC11_IC_val, G, GV, US, CS)

  if (.not.restart .or. (CS%tracers_may_reinit .and. &
      .not.query_initialized(CS%CFC12, CS%CFC12_name, CS%restart_CSp))) &
    call init_tracer_CFC(h, CS%CFC12, CS%CFC12_name, CS%CFC12_land_val, &
                         CS%CFC12_IC_val, G, GV, US, CS)


  ! cmor diagnostics
  ! CFC11 cmor conventions: http://clipc-services.ceda.ac.uk/dreq/u/42625c97b8fe75124a345962c4430982.html
  CS%id_cfc11_cmor = register_diag_field('ocean_model', 'cfc11', diag%axesTL, day,   &
    'Mole Concentration of CFC11 in Sea Water', 'mol m-3')
  ! CFC12 cmor conventions: http://clipc-services.ceda.ac.uk/dreq/u/3ab8e10027d7014f18f9391890369235.html
  CS%id_cfc12_cmor = register_diag_field('ocean_model', 'cfc12', diag%axesTL, day,   &
    'Mole Concentration of CFC12 in Sea Water', 'mol m-3')

  if (associated(OBC)) then
  ! Steal from updated DOME in the fullness of time.
  ! GMM: TODO this must be coded if we intend to use this module in regional applications
  endif

end subroutine initialize_CFC_cap

!>This subroutine initializes a tracer array.
subroutine init_tracer_CFC(h, tr, name, land_val, IC_val, G, GV, US, CS)
  type(ocean_grid_type),                     intent(in)  :: G        !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV       !< The ocean's vertical grid structure.
  type(unit_scale_type),                     intent(in)  :: US       !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h        !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: tr       !< The tracer concentration array
  character(len=*),                          intent(in)  :: name     !< The tracer name
  real,                                      intent(in)  :: land_val !< A value the tracer takes over land
  real,                                      intent(in)  :: IC_val   !< The initial condition value for the tracer
  type(CFC_cap_CS),                          pointer     :: CS       !< The control structure returned by a
                                                                     !! previous call to register_CFC_cap.

  ! local variables
  logical :: OK
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (len_trim(CS%IC_file) > 0) then
    !  Read the tracer concentrations from a netcdf file.
    if (.not.file_exists(CS%IC_file, G%Domain)) &
      call MOM_error(FATAL, "initialize_CFC_cap: Unable to open "//CS%IC_file)
    if (CS%Z_IC_file) then
      OK = tracer_Z_init(tr, h, CS%IC_file, name, G, GV, US)
      if (.not.OK) then
        OK = tracer_Z_init(tr, h, CS%IC_file, trim(name), G, GV, US)
        if (.not.OK) call MOM_error(FATAL,"initialize_CFC_cap: "//&
                "Unable to read "//trim(name)//" from "//&
                trim(CS%IC_file)//".")
      endif
    else
      call MOM_read_data(CS%IC_file, trim(name), tr, G%Domain)
    endif
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      if (G%mask2dT(i,j) < 0.5) then
        tr(i,j,k) = land_val
      else
        tr(i,j,k) = IC_val
      endif
    enddo ; enddo ; enddo
  endif

end subroutine init_tracer_CFC

!> Applies diapycnal diffusion, souces and sinks and any other column
!! tracer physics to the CFC cap tracers. CFCs are relatively simple,
!! as they are passive tracers with only a surface flux as a source.
subroutine CFC_cap_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G     !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV    !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: ea    !< an array to which the amount of fluid entrained
                                               !! from the layer above during this call will be
                                               !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: eb    !< an array to which the amount of fluid entrained
                                               !! from the layer below during this call will be
                                               !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes!< A structure containing pointers to thermodynamic
                                               !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt    !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US    !< A dimensional unit scaling type
  type(CFC_cap_CS),        pointer    :: CS    !< The control structure returned by a
                                               !! previous call to register_CFC_cap.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                               !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                               !! fluxes can be applied [H ~> m or kg m-2]

  ! The arguments to this subroutine are redundant in that
  !     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  integer :: i, j, k, m, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return

  ! Use a tridiagonal solver to determine the concentrations after the
  ! surface source is applied and diapycnal advection and diffusion occurs.
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do k=1,nz ;do j=js,je ; do i=is,ie
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo
    call applyTracerBoundaryFluxesInOut(G, GV, CS%CFC11, dt, fluxes, h_work, &
                                        evap_CFL_limit, minimum_forcing_depth)
    call tracer_vertdiff(h_work, ea, eb, dt, CS%CFC11, G, GV, sfc_flux=fluxes%cfc11_flux)

    do k=1,nz ;do j=js,je ; do i=is,ie
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo
    call applyTracerBoundaryFluxesInOut(G, GV, CS%CFC12, dt, fluxes, h_work, &
                                        evap_CFL_limit, minimum_forcing_depth)
    call tracer_vertdiff(h_work, ea, eb, dt, CS%CFC12, G, GV, sfc_flux=fluxes%cfc12_flux)
  else
    call tracer_vertdiff(h_old, ea, eb, dt, CS%CFC11, G, GV, sfc_flux=fluxes%cfc11_flux)
    call tracer_vertdiff(h_old, ea, eb, dt, CS%CFC12, G, GV, sfc_flux=fluxes%cfc12_flux)
  endif

  ! If needed, write out any desired diagnostics from tracer sources & sinks here.
  if (CS%id_cfc11_cmor > 0) call post_data(CS%id_cfc11_cmor, (GV%Rho0*US%R_to_kg_m3)*CS%CFC11, CS%diag)
  if (CS%id_cfc12_cmor > 0) call post_data(CS%id_cfc12_cmor, (GV%Rho0*US%R_to_kg_m3)*CS%CFC12, CS%diag)

end subroutine CFC_cap_column_physics

!> Calculates the mass-weighted integral of all tracer stocks,
!! returning the number of stocks it has calculated.  If the stock_index
!! is present, only the stock corresponding to that coded index is returned.
function CFC_cap_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),           intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type),         intent(in)    :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                   intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2].
  type(EFP_type), dimension(:),    intent(out)   :: stocks !< The mass-weighted integrated amount of each
                                                           !! tracer, in kg times concentration units [kg conc]
  type(CFC_cap_CS),                pointer       :: CS     !< The control structure returned by a
                                                           !! previous call to register_CFC_cap.
  character(len=*), dimension(:),  intent(out)   :: names  !< The names of the stocks calculated.
  character(len=*), dimension(:),  intent(out)   :: units  !< The units of the stocks calculated.
  integer, optional,               intent(in)    :: stock_index !< The coded index of a specific
                                                                !! stock being sought.
  integer                                        :: CFC_cap_stock !< The number of stocks calculated here.


  CFC_cap_stock = 0
  if (.not.associated(CS)) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  call query_vardesc(CS%CFC11_desc, name=names(1), units=units(1), caller="CFC_cap_stock")
  call query_vardesc(CS%CFC12_desc, name=names(2), units=units(2), caller="CFC_cap_stock")
  units(1) = trim(units(1))//" kg" ; units(2) = trim(units(2))//" kg"

  stocks(1) = global_mass_int_EFP(h, G, GV, CS%CFC11, on_PE_only=.true.)
  stocks(2) = global_mass_int_EFP(h, G, GV, CS%CFC12, on_PE_only=.true.)

  CFC_cap_stock = 2

end function CFC_cap_stock

!> Extracts the ocean surface CFC concentrations and copies them to sfc_state.
subroutine CFC_cap_surface_state(sfc_state, G, CS)
  type(ocean_grid_type),   intent(in)    :: G !< The ocean's grid structure.
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                              !! describe the surface state of the ocean.
  type(CFC_cap_CS),        pointer       :: CS!< The control structure returned by a previous
                                              !! call to register_CFC_cap.

  ! Local variables
  integer :: i, j, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.not.associated(CS)) return

  do j=js,je ; do i=is,ie
    sfc_state%sfc_cfc11(i,j) = CS%CFC11(i,j,1)
    sfc_state%sfc_cfc12(i,j) = CS%CFC12(i,j,1)
  enddo ; enddo

end subroutine CFC_cap_surface_state

!> Orchestrates the calculation of the CFC fluxes [mol m-2 s-1], including getting the ATM
!! concentration, and calculating the solubility, Schmidt number, and gas exchange.
subroutine CFC_cap_fluxes(fluxes, sfc_state, G, US, Rho0, Time, id_cfc11_atm, id_cfc12_atm)
  type(ocean_grid_type),        intent(in   ) :: G  !< The ocean's grid structure.
  type(unit_scale_type),        intent(in  )  :: US !< A dimensional unit scaling type
  type(surface),                intent(in   ) :: sfc_state !< A structure containing fields
                                              !! that describe the surface state of the ocean.
  type(forcing),                intent(inout) :: fluxes !< A structure containing pointers
                                              !! to thermodynamic and tracer forcing fields. Unused fields
                                              !! have NULL ptrs.
  real,                         intent(in   ) :: Rho0 !< The mean ocean density [R ~> kg m-3]
  type(time_type),              intent(in   ) :: Time !< The time of the fluxes, used for interpolating the
                                              !! CFC's concentration in the atmosphere.
  integer,           optional,  intent(inout):: id_cfc11_atm !< id number for time_interp_external.
  integer,           optional,  intent(inout):: id_cfc12_atm !< id number for time_interp_external.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    kw_wo_sc_no_term, &  ! gas transfer velocity, without the Schmidt number term [Z T-1 ~> m s-1].
    kw, &           ! gas transfer velocity [Z T-1 ~> m s-1].
    cair, &         ! The surface gas concentration in equilibrium with the atmosphere (saturation concentration)
                    ! [mol kg-1].
    cfc11_atm,     & !< CFC11 concentration in the atmopshere [pico mol/mol]
    cfc12_atm        !< CFC11 concentration in the atmopshere [pico mol/mol]
  real :: ta        ! Absolute sea surface temperature [hectoKelvin]
  real :: sal       ! Surface salinity [PSU].
  real :: alpha_11  ! The solubility of CFC 11 [mol kg-1 atm-1].
  real :: alpha_12  ! The solubility of CFC 12 [mol kg-1 atm-1].
  real :: sc_11, sc_12 ! The Schmidt numbers of CFC 11 and CFC 12 [nondim].
  real :: kw_coeff     ! A coefficient used to compute the piston velocity [Z T-1 T2 L-2 = Z T L-2 ~> s / m]
  real, parameter :: pa_to_atm = 9.8692316931427e-6 ! factor for converting from Pa to atm [atm Pa-1].
  real :: press_to_atm ! converts from model pressure units to atm [atm T2 R-1 L-2 ~> atm Pa-1]
  integer :: i, j, m, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  ! CFC11 ATM concentration
  if (present(id_cfc11_atm) .and. (id_cfc11_atm /= -1)) then
    call time_interp_external(id_cfc11_atm, Time, cfc11_atm)
    ! convert from ppt (pico mol/mol) to mol/mol
    cfc11_atm = cfc11_atm * 1.0e-12
  else
    ! TODO: create cfc11_atm internally
    call MOM_error(FATAL, "CFC_cap_fluxes: option to create cfc11_atm internally" //&
                          "has not been implemented yet.")
  endif

  ! CFC12 ATM concentration
  if (present(id_cfc12_atm) .and. (id_cfc12_atm /= -1)) then
    call time_interp_external(id_cfc12_atm, Time, cfc12_atm)
    ! convert from ppt (pico mol/mol) to mol/mol
    cfc12_atm = cfc12_atm * 1.0e-12
  else
    ! TODO: create cfc11_atm internally
    call MOM_error(FATAL, "CFC_cap_fluxes: option to create cfc12_atm internally" //&
                          "has not been implemented yet.")
  endif

  !---------------------------------------------------------------------
  !     Gas exchange/piston velocity parameter
  !---------------------------------------------------------------------
  ! From a = 0.251 cm/hr s^2/m^2 in Wannikhof 2014
  !        = 6.97e-7 m/s s^2/m^2 [Z T-1 T2 L-2 = Z T L-2 ~> s / m]
  kw_coeff = (US%m_to_Z*US%s_to_T*US%L_to_m**2) * 6.97e-7

  ! set unit conversion factors
  press_to_atm = US%R_to_kg_m3*US%L_T_to_m_s**2 * pa_to_atm

  do j=js,je ; do i=is,ie
    ! ta in hectoKelvin
    ta = max(0.01, (sfc_state%SST(i,j) + 273.15) * 0.01)
    sal = sfc_state%SSS(i,j)

    ! Calculate solubilities
    call get_solubility(alpha_11, alpha_12, ta, sal , G%mask2dT(i,j))

    ! Calculate Schmidt numbers using coefficients given by
    ! Wanninkhof (2014); doi:10.4319/lom.2014.12.351.
    call comp_CFC_schmidt(sfc_state%SST(i,j), sc_11, sc_12)

    kw_wo_sc_no_term(i,j) = kw_coeff * ((1.0 - fluxes%ice_fraction(i,j))*fluxes%u10_sqr(i,j))

    ! air concentrations and cfcs BC's fluxes
    ! CFC flux units: CU R Z T-1 = mol kg-1 R Z T-1 ~> mol m-2 s-1
    kw(i,j) = kw_wo_sc_no_term(i,j) * sqrt(660.0 / sc_11)
    cair(i,j) = press_to_atm * alpha_11 * cfc11_atm(i,j) * fluxes%p_surf_full(i,j)
    fluxes%cfc11_flux(i,j) = kw(i,j) * (cair(i,j) - sfc_state%sfc_CFC11(i,j)) * Rho0

    kw(i,j) = kw_wo_sc_no_term(i,j) * sqrt(660.0 / sc_12)
    cair(i,j) = press_to_atm * alpha_12 * cfc12_atm(i,j) * fluxes%p_surf_full(i,j)
    fluxes%cfc12_flux(i,j) = kw(i,j) * (cair(i,j) - sfc_state%sfc_CFC12(i,j)) * Rho0
  enddo ; enddo

end subroutine CFC_cap_fluxes

!> Calculates the CFC's solubility function following Warner and Weiss (1985) DSR, vol 32.
subroutine get_solubility(alpha_11, alpha_12, ta, sal , mask)
  real, intent(inout) :: alpha_11 !< The solubility of CFC 11 [mol kg-1 atm-1]
  real, intent(inout) :: alpha_12 !< The solubility of CFC 12 [mol kg-1 atm-1]
  real, intent(in   ) :: ta       !< Absolute sea surface temperature [hectoKelvin]
  real, intent(in   ) :: sal      !< Surface salinity [PSU].
  real, intent(in   ) :: mask     !< ocean mask [nondim]

  ! Local variables

  ! Coefficients for calculating CFC11 solubilities
  ! from Table 5 in Warner and Weiss (1985) DSR, vol 32.

  real, parameter :: d1_11 = -232.0411    ! [nondim]
  real, parameter :: d2_11 =  322.5546    ! [hectoKelvin-1]
  real, parameter :: d3_11 =  120.4956    ! [log(hectoKelvin)-1]
  real, parameter :: d4_11 =   -1.39165   ! [hectoKelvin-2]

  real, parameter :: e1_11 =   -0.146531  ! [PSU-1]
  real, parameter :: e2_11 =    0.093621  ! [PSU-1 hectoKelvin-1]
  real, parameter :: e3_11 =   -0.0160693 ! [PSU-2 hectoKelvin-2]

  ! Coefficients for calculating CFC12 solubilities
  ! from Table 5 in Warner and Weiss (1985) DSR, vol 32.

  real, parameter :: d1_12 = -220.2120    ! [nondim]
  real, parameter :: d2_12 =  301.8695    ! [hectoKelvin-1]
  real, parameter :: d3_12 =  114.8533    ! [log(hectoKelvin)-1]
  real, parameter :: d4_12 =   -1.39165   ! [hectoKelvin-2]

  real, parameter :: e1_12 =   -0.147718  ! [PSU-1]
  real, parameter :: e2_12 =    0.093175  ! [PSU-1 hectoKelvin-1]
  real, parameter :: e3_12 =   -0.0157340 ! [PSU-2 hectoKelvin-2]

  real :: factor ! introduce units to result [mol kg-1 atm-1]

  ! Eq. 9 from Warner and Weiss (1985) DSR, vol 32.
  factor = 1.0
  alpha_11 = exp(d1_11 + d2_11/ta + d3_11*log(ta) + d4_11*ta**2 +&
                 sal * ((e3_11 * ta + e2_11) * ta + e1_11)) * &
             factor * mask
  alpha_12 = exp(d1_12 + d2_12/ta + d3_12*log(ta) + d4_12*ta**2 +&
                 sal * ((e3_12 * ta + e2_12) * ta + e1_12)) * &
             factor * mask

end subroutine get_solubility


!> Compute Schmidt numbers of CFCs following Wanninkhof (2014); doi:10.4319/lom.2014.12.351
!! Range of validity of fit is -2:40.
subroutine comp_CFC_schmidt(sst_in, cfc11_sc, cfc12_sc)
  real, intent(in)    :: sst_in   !< The sea surface temperature [degC].
  real, intent(inout) :: cfc11_sc !< Schmidt number of CFC11 [nondim].
  real, intent(inout) :: cfc12_sc !< Schmidt number of CFC12 [nondim].

  !local variables
  real , parameter :: a_11 = 3579.2    ! CFC11 Schmidt number fit coefficient [nondim]
  real , parameter :: b_11 = -222.63   ! CFC11 Schmidt number fit coefficient [degC-1]
  real , parameter :: c_11 = 7.5749    ! CFC11 Schmidt number fit coefficient [degC-2]
  real , parameter :: d_11 = -0.14595  ! CFC11 Schmidt number fit coefficient [degC-3]
  real , parameter :: e_11 = 0.0011874 ! CFC11 Schmidt number fit coefficient [degC-4]
  real , parameter :: a_12 = 3828.1    ! CFC12 Schmidt number fit coefficient [nondim]
  real , parameter :: b_12 = -249.86   ! CFC12 Schmidt number fit coefficient [degC-1]
  real , parameter :: c_12 = 8.7603    ! CFC12 Schmidt number fit coefficient [degC-2]
  real , parameter :: d_12 = -0.1716   ! CFC12 Schmidt number fit coefficient [degC-3]
  real , parameter :: e_12 = 0.001408  ! CFC12 Schmidt number fit coefficient [degC-4]
  real             :: sst  ! A range-limited sea surface temperature [degC]


  ! clip SST to avoid bad values
  sst = MAX(-2.0, MIN(40.0, sst_in))
  cfc11_sc = a_11 + sst * (b_11 + sst * (c_11 + sst * (d_11 + sst * e_11)))
  cfc12_sc = a_12 + sst * (b_12 + sst * (c_12 + sst * (d_12 + sst * e_12)))

end subroutine comp_CFC_schmidt

!> Deallocate any memory associated with the CFC cap tracer package
subroutine CFC_cap_end(CS)
  type(CFC_cap_CS), pointer :: CS !< The control structure returned by a
                                  !! previous call to register_CFC_cap.
  ! local variables
  integer :: m

  if (associated(CS)) then
    if (associated(CS%CFC11)) deallocate(CS%CFC11)
    if (associated(CS%CFC12)) deallocate(CS%CFC12)

    deallocate(CS)
  endif
end subroutine CFC_cap_end

!> Unit tests for the CFC cap module.
logical function CFC_cap_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, output additional
                                 !! information for debugging unit tests

  ! Local variables
  real               :: dummy1, dummy2, ta, sal
  character(len=120) :: test_name ! Title of the unit test

  CFC_cap_unit_tests = .false.
  write(stdout,*) '==== MOM_CFC_cap ======================='

  ! test comp_CFC_schmidt, Table 1 in Wanninkhof (2014); doi:10.4319/lom.2014.12.351
  test_name = 'Schmidt number calculation'
  call comp_CFC_schmidt(20.0, dummy1, dummy2)
  CFC_cap_unit_tests = CFC_cap_unit_tests .or. &
                       compare_values(verbose, test_name, dummy1, 1179.0, 0.5)
  CFC_cap_unit_tests = CFC_cap_unit_tests .or. &
                       compare_values(verbose, test_name, dummy2, 1188.0, 0.5)

  if (.not. CFC_cap_unit_tests) write(stdout,'(2x,a)') "Passed "//test_name

  test_name = 'Solubility function, SST = 1.0 C, and SSS = 10 psu'
  ta = max(0.01, (1.0 + 273.15) * 0.01); sal = 10.
  ! cfc1 = 3.238 10-2 mol kg-1 atm-1
  ! cfc2 = 7.943 10-3 mol kg-1 atm-1
  call get_solubility(dummy1, dummy2, ta, sal , 1.0)
  CFC_cap_unit_tests = CFC_cap_unit_tests .or. &
                       compare_values(verbose, test_name, dummy1, 3.238e-2, 5.0e-6)
  CFC_cap_unit_tests = CFC_cap_unit_tests .or. &
                       compare_values(verbose, test_name, dummy2, 7.943e-3, 5.0e-6)

  if (.not. CFC_cap_unit_tests) write(stdout,'(2x,a)')"Passed "//test_name

  test_name = 'Solubility function, SST = 20.0 C, and SSS = 35 psu'
  ta = max(0.01, (20.0 + 273.15) * 0.01); sal = 35.
  ! cfc1 = 0.881 10-2 mol kg-1 atm-1
  ! cfc2 = 2.446 10-3 mol kg-1 atm-1
  call get_solubility(dummy1, dummy2, ta, sal , 1.0)
  CFC_cap_unit_tests = CFC_cap_unit_tests .or. &
                       compare_values(verbose, test_name, dummy1, 8.8145e-3, 5.0e-8)
  CFC_cap_unit_tests = CFC_cap_unit_tests .or. &
                       compare_values(verbose, test_name, dummy2, 2.4462e-3, 5.0e-8)
  if (.not. CFC_cap_unit_tests) write(stdout,'(2x,a)')"Passed "//test_name

end function CFC_cap_unit_tests

!> Test that ans and calc are approximately equal by computing the difference
!! and comparing it against limit.
logical function compare_values(verbose, test_name, calc, ans, limit)
  logical,             intent(in) :: verbose   !< If true, write results to stdout
  character(len=80),   intent(in) :: test_name !< Brief description of the unit test
  real,                intent(in) :: calc      !< computed value
  real,                intent(in) :: ans       !< correct value
  real,                intent(in) :: limit     !< value above which test fails

  ! Local variables
  real :: diff

  diff = ans - calc

  compare_values = .false.
  if (diff > limit ) then
    compare_values = .true.
    write(stdout,*) "CFC_cap_unit_tests, UNIT TEST FAILED: ", test_name
    write(stdout,10) calc, ans
  elseif (verbose) then
    write(stdout,10) calc, ans
  endif

10 format("calc=",f20.16," ans",f20.16)
end function compare_values

!> \namespace mom_CFC_cap
!!
!!     This module contains the code that is needed to simulate
!!   CFC-11 and CFC-12 using atmospheric and sea ice variables
!!   provided via cap (only NUOPC cap is implemented so far).

end module MOM_CFC_cap
