!> A sample tracer package that has striped initial conditions
module USER_tracer_example

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coupler_types, only : set_coupler_type_data, atmos_ocn_coupler_flux
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_io, only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public USER_register_tracer_example, USER_initialize_tracer, USER_tracer_stock
public tracer_column_physics, USER_tracer_surface_state, USER_tracer_example_end

integer, parameter :: NTR = 1 !< The number of tracers in this module.

!> The control structure for the USER_tracer_example module
type, public :: USER_tracer_example_CS ; private
  logical :: coupled_tracers = .false. !< These tracers are not offered to the coupler.
  character(len=200) :: tracer_IC_file !< The full path to the IC file, or " "
                                       !! to initialize internally.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL()  !< The array of tracers used in this subroutine, in g m-3?
  real :: land_val(NTR) = -1.0 !< The value of tr that is used where land is masked out.
  logical :: use_sponge    !< If true, sponges may be applied somewhere in the domain.

  integer, dimension(NTR) :: ind_tr !< Indices returned by atmos_ocn_coupler_flux if it is used and the
                                    !! surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the timing of diagnostic output.

  type(vardesc) :: tr_desc(NTR) !< Descriptions of each of the tracers.
end type USER_tracer_example_CS

contains

!> This subroutine is used to register tracer fields and subroutines to be used with MOM.
function USER_register_tracer_example(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),    intent(in)   :: HI   !< A horizontal index type structure
  type(verticalGrid_type), intent(in)   :: GV   !< The ocean's vertical grid structure
  type(param_file_type),   intent(in)   :: param_file !< A structure to parse for run-time parameters
  type(USER_tracer_example_CS), pointer :: CS   !< A pointer that is set to point to the control
                                                !! structure for this module
  type(tracer_registry_type), pointer   :: tr_Reg !< A pointer that is set to point to the control
                                                  !! structure for the tracer advection and
                                                  !! diffusion module
  type(MOM_restart_CS),       pointer   :: restart_CS !< A pointer to the restart control structure

! Local variables
  character(len=80)  :: name, longname
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "tracer_example" ! This module's name.
  character(len=200) :: inputdir
  character(len=48) :: flux_units ! The units for tracer fluxes, usually
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: USER_register_tracer_example
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "USER_register_tracer_example called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "TRACER_EXAMPLE_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial "//&
                 "conditions for the DOME tracers, or blank to initialize "//&
                 "them internally.", default=" ")
  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    CS%tracer_IC_file = trim(slasher(inputdir))//trim(CS%tracer_IC_file)
    call log_param(param_file, mdl, "INPUTDIR/TRACER_EXAMPLE_IC_FILE", &
                   CS%tracer_IC_file)
  endif
  call get_param(param_file, mdl, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified from MOM_initialization.F90.", default=.false.)

  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR)) ; CS%tr(:,:,:,:) = 0.0

  do m=1,NTR
    if (m < 10) then ; write(name,'("tr",I1.1)') m
    else ; write(name,'("tr",I2.2)') m ; endif
    write(longname,'("Concentration of Tracer ",I2.2)') m
    CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)

    ! This needs to be changed if the units of tracer are changed above.
    if (GV%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
    else ; flux_units = "kg s-1" ; endif

    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                         name=name, longname=longname, units="kg kg-1", &
                         registry_diags=.true., flux_units=flux_units, &
                         restart_CS=restart_CS)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = atmos_ocn_coupler_flux(trim(name)//'_flux', &
          flux_type=' ', implementation=' ', caller="USER_register_tracer_example")
  enddo

  CS%tr_Reg => tr_Reg
  USER_register_tracer_example = .true.
end function USER_register_tracer_example

!> This subroutine initializes the NTR tracer fields in tr(:,:,:,:)
!! and it sets up the tracer output.
subroutine USER_initialize_tracer(restart, day, G, GV, h, diag, OBC, CS, &
                                  sponge_CSp)
  logical,                            intent(in) :: restart !< .true. if the fields have already
                                                         !! been read from a restart file.
  type(time_type),            target, intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                      intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl),            target, intent(in) :: diag !< A structure that is used to regulate
                                                         !! diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(USER_tracer_example_CS),       pointer    :: CS   !< The control structure returned by a previous
                                                         !! call to USER_register_tracer_example.
  type(sponge_CS),                    pointer    :: sponge_CSp    !< A pointer to the control structure
                                                                  !! for the sponges, if they are in use.

! Local variables
  real, allocatable :: temp(:,:,:)
  character(len=32) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for tracer fluxes, usually
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  real, pointer :: tr_ptr(:,:,:) => NULL()
  real :: PI     ! 3.1415926... calculated as 4*atan(1)
  real :: tr_y   ! Initial zonally uniform tracer concentrations.
  real :: dist2  ! The distance squared from a line [m2].
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB, lntr

  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  lntr = NTR ! Avoids compile-time warning when NTR<2
  CS%Time => day
  CS%diag => diag

  if (.not.restart) then
    if (len_trim(CS%tracer_IC_file) >= 1) then
!  Read the tracer concentrations from a netcdf file.
      if (.not.file_exists(CS%tracer_IC_file, G%Domain)) &
        call MOM_error(FATAL, "USER_initialize_tracer: Unable to open "// &
                        CS%tracer_IC_file)
      do m=1,NTR
        call query_vardesc(CS%tr_desc(m), name, caller="USER_initialize_tracer")
        call MOM_read_data(CS%tracer_IC_file, trim(name), CS%tr(:,:,:,m), G%Domain)
      enddo
    else
      do m=1,NTR
        do k=1,nz ; do j=js,je ; do i=is,ie
          CS%tr(i,j,k,m) = 1.0e-20 ! This could just as well be 0.
        enddo ; enddo ; enddo
      enddo

!    This sets a stripe of tracer across the basin.
      PI = 4.0*atan(1.0)
      do j=js,je
        dist2 = (G%Rad_Earth * PI / 180.0)**2 * &
               (G%geoLatT(i,j) - 40.0) * (G%geoLatT(i,j) - 40.0)
        tr_y = 0.5*exp(-dist2/(1.0e5*1.0e5))

        do k=1,nz ; do i=is,ie
!      This adds the stripes of tracer to every layer.
          CS%tr(i,j,k,1) = CS%tr(i,j,k,1) + tr_y
        enddo ; enddo
      enddo
    endif
  endif ! restart

  if ( CS%use_sponge ) then
!   If sponges are used, this example damps tracers in sponges in the
! northern half of the domain to 1 and tracers in the southern half
! to 0.  For any tracers that are not damped in the sponge, the call
! to set_up_sponge_field can simply be omitted.
    if (.not.associated(sponge_CSp)) &
      call MOM_error(FATAL, "USER_initialize_tracer: "// &
        "The pointer to sponge_CSp must be associated if SPONGE is defined.")

    allocate(temp(G%isd:G%ied,G%jsd:G%jed,nz))
    do k=1,nz ; do j=js,je ; do i=is,ie
      if (G%geoLatT(i,j) > 700.0 .and. (k > nz/2)) then
        temp(i,j,k) = 1.0
      else
        temp(i,j,k) = 0.0
      endif
    enddo ; enddo ; enddo

!   do m=1,NTR
    do m=1,1
      ! This is needed to force the compiler not to do a copy in the sponge
      ! calls.  Curses on the designers and implementers of Fortran90.
      tr_ptr => CS%tr(:,:,:,m)
      call set_up_sponge_field(temp, tr_ptr, G, GV, nz, sponge_CSp)
    enddo
    deallocate(temp)
  endif

  if (associated(OBC)) then
    call query_vardesc(CS%tr_desc(1), name, caller="USER_initialize_tracer")
    if (OBC%specified_v_BCs_exist_globally) then
      ! Steal from updated DOME in the fullness of time.
    else
      ! Steal from updated DOME in the fullness of time.
    endif
    ! All tracers but the first have 0 concentration in their inflows. As this
    ! is the default value, the following calls are unnecessary.
    do m=2,lntr
      call query_vardesc(CS%tr_desc(m), name, caller="USER_initialize_tracer")
      ! Steal from updated DOME in the fullness of time.
    enddo
  endif

end subroutine USER_initialize_tracer

!> This subroutine applies diapycnal diffusion and any other column
!! tracer physics or chemistry to the tracers from this file.
!! This is a simple example of a set of advected passive tracers.
!! The arguments to this subroutine are redundant in that
!!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)
subroutine tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(USER_tracer_example_CS), pointer :: CS !< The control structure returned by a previous
                                              !! call to USER_register_tracer_example.

! Local variables
  real :: hold0(SZI_(G))       ! The original topmost layer thickness,
                               ! with surface mass fluxes added back, m.
  real :: b1(SZI_(G))          ! b1 and c1 are variables used by the
  real :: c1(SZI_(G),SZK_(GV)) ! tridiagonal solver.
  real :: d1(SZI_(G))          ! d1=1-c1 is used by the tridiagonal solver.
  real :: h_neglect            ! A thickness that is so small it is usually lost
                               ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: b_denom_1 ! The first term in the denominator of b1 [H ~> m or kg m-2].
  integer :: i, j, k, is, ie, js, je, nz, m

! The following array (trdc) determines the behavior of the tracer
! diapycnal advection.  The first element is 1 if tracers are
! passively advected.  The second and third are the concentrations
! to which downwelling and upwelling water are set, respectively.
! For most (normal) tracers, the appropriate vales are {1,0,0}.

  real :: trdc(3)
!   Uncomment the following line to dye both upwelling and downwelling.
! data trdc / 0.0,1.0,1.0 /
!   Uncomment the following line to dye downwelling.
! data trdc / 0.0,1.0,0.0 /
!   Uncomment the following line to dye upwelling.
! data trdc / 0.0,0.0,1.0 /
!   Uncomment the following line for tracer concentrations to be set
! to zero in any diapycnal motions.
! data trdc / 0.0,0.0,0.0 /
!   Uncomment the following line for most "physical" tracers, which
! are advected diapycnally in the usual manner.
  data trdc / 1.0,0.0,0.0 /
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  h_neglect = GV%H_subroundoff

  do j=js,je
    do i=is,ie
!   The following line is appropriate for quantities like salinity
! that are left behind by evaporation, and any surface fluxes would
! be explicitly included in the flux structure.
      hold0(i) = h_old(i,j,1)
!   The following line is appropriate for quantities like temperature
! that can be assumed to have the same concentration in evaporation
! as they had in the water.  The explicit surface fluxes here would
! reflect differences in concentration from the ambient water, not
! the absolute fluxes.
  !   hold0(i) = h_old(i,j,1) + ea(i,j,1)
      b_denom_1 = h_old(i,j,1) + ea(i,j,1) + h_neglect
      b1(i) = 1.0 / (b_denom_1 + eb(i,j,1))
!       d1(i) = b_denom_1 * b1(i)
      d1(i) = trdc(1) * (b_denom_1 * b1(i)) + (1.0 - trdc(1))
      do m=1,NTR
        CS%tr(i,j,1,m) = b1(i)*(hold0(i)*CS%tr(i,j,1,m) + trdc(3)*eb(i,j,1))
 !      Add any surface tracer fluxes to the preceding line.
      enddo
    enddo
    do k=2,nz ; do i=is,ie
      c1(i,k) = trdc(1) * eb(i,j,k-1) * b1(i)
      b_denom_1 = h_old(i,j,k) + d1(i)*ea(i,j,k) + h_neglect
      b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
      d1(i) = trdc(1) * (b_denom_1 * b1(i)) + (1.0 - trdc(1))
      do m=1,NTR
        CS%tr(i,j,k,m) = b1(i) * (h_old(i,j,k)*CS%tr(i,j,k,m) + &
                 ea(i,j,k)*(trdc(1)*CS%tr(i,j,k-1,m)+trdc(2)) + &
                 eb(i,j,k)*trdc(3))
      enddo
    enddo ; enddo
    do m=1,NTR ; do k=nz-1,1,-1 ; do i=is,ie
      CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + c1(i,k+1)*CS%tr(i,j,k+1,m)
    enddo ; enddo ; enddo
  enddo

end subroutine tracer_column_physics

!> This function calculates the mass-weighted integral of all tracer stocks,
!! returning the number of stocks it has calculated.  If the stock_index
!! is present, only the stock corresponding to that coded index is returned.
function USER_tracer_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                      intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(:),                 intent(out)   :: stocks !< the mass-weighted integrated amount of each
                                                              !! tracer, in kg times concentration units [kg conc].
  type(USER_tracer_example_CS),       pointer       :: CS     !< The control structure returned by a
                                                              !! previous call to register_USER_tracer.
  character(len=*), dimension(:),     intent(out)   :: names  !< The names of the stocks calculated.
  character(len=*), dimension(:),     intent(out)   :: units  !< The units of the stocks calculated.
  integer, optional,                  intent(in)    :: stock_index !< The coded index of a specific stock
                                                              !! being sought.
  integer                                           :: USER_tracer_stock !< Return value: the number of
                                                              !! stocks calculated here.

  ! Local variables
  real :: stock_scale ! The dimensional scaling factor to convert stocks to kg [kg H-1 L-2 ~> kg m-3 or nondim]
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  USER_tracer_stock = 0
  if (.not.associated(CS)) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  stock_scale = G%US%L_to_m**2 * GV%H_to_kg_m2
  do m=1,NTR
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="USER_tracer_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      stocks(m) = stocks(m) + CS%tr(i,j,k,m) * (G%mask2dT(i,j) * G%areaT(i,j) * h(i,j,k))
    enddo ; enddo ; enddo
    stocks(m) = stock_scale * stocks(m)
  enddo
  USER_tracer_stock = NTR

end function USER_tracer_stock

!> This subroutine extracts the surface fields from this tracer package that
!! are to be shared with the atmosphere in coupled configurations.
subroutine USER_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),        intent(in)    :: G     !< The ocean's grid structure
  type(verticalGrid_type),      intent(in)    :: GV    !< The ocean's vertical grid structure
  type(surface),                intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(USER_tracer_example_CS), pointer       :: CS !< The control structure returned by a previous
                                                    !! call to register_USER_tracer.

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

  integer :: m, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,ntr
      !   This call loads the surface values into the appropriate array in the
      ! coupler-type structure.
      call set_coupler_type_data(CS%tr(:,:,1,m), CS%ind_tr(m), sfc_state%tr_fields, &
                   idim=(/isd, is, ie, ied/), jdim=(/jsd, js, je, jed/) )
    enddo
  endif

end subroutine USER_tracer_surface_state

!> Clean up allocated memory at the end.
subroutine USER_tracer_example_end(CS)
  type(USER_tracer_example_CS), pointer :: CS !< The control structure returned by a previous
                                              !! call to register_USER_tracer.
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end subroutine USER_tracer_example_end

!> \namespace user_tracer_example
!!
!!  Original by Robert Hallberg, 2002
!!
!!    This file contains an example of the code that is needed to set
!!  up and use a set (in this case one) of dynamically passive tracers.
!!
!!    A single subroutine is called from within each file to register
!!  each of the tracers for reinitialization and advection and to
!!  register the subroutine that initializes the tracers and set up
!!  their output and the subroutine that does any tracer physics or
!!  chemistry along with diapycnal mixing (included here because some
!!  tracers may float or swim vertically or dye diapycnal processes).
end module USER_tracer_example
