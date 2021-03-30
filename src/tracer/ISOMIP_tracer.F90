!> Routines used to set up and use a set of (one for now)
!! dynamically passive tracers in the ISOMIP configuration.
!!
!! For now, just one passive tracer is injected in
!! the sponge layer.
module ISOMIP_tracer

! This file is part of MOM6. See LICENSE.md for the license.

! Original sample tracer package by Robert Hallberg, 2002
! Adapted to the ISOMIP test case by Gustavo Marques, May 2016

use MOM_coms, only : max_across_PEs
use MOM_coupler_types, only : set_coupler_type_data, atmos_ocn_coupler_flux
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_hor_index, only : hor_index_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : MOM_restart_CS
use MOM_ALE_sponge, only : set_up_ALE_sponge_field, ALE_sponge_CS
use MOM_time_manager, only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!< Publicly available functions
public register_ISOMIP_tracer, initialize_ISOMIP_tracer
public ISOMIP_tracer_column_physics, ISOMIP_tracer_surface_state, ISOMIP_tracer_end

integer, parameter :: ntr = 1 !< ntr is the number of tracers in this module.

!> ISOMIP tracer package control structure
type, public :: ISOMIP_tracer_CS ; private
  logical :: coupled_tracers = .false.  !< These tracers are not offered to the coupler.
  character(len = 200) :: tracer_IC_file !< The full path to the IC file, or " " to initialize internally.
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the MOM tracer registry
  real, pointer :: tr(:,:,:,:) => NULL()   !< The array of tracers used in this package, in g m-3?
  real :: land_val(NTR) = -1.0 !< The value of tr used where land is masked out.
  logical :: use_sponge    !< If true, sponges may be applied somewhere in the domain.

  integer, dimension(NTR) :: ind_tr !< Indices returned by atmos_ocn_coupler_flux
             !< if it is used and the surface tracer concentrations are to be
             !< provided to the coupler.

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  type(vardesc) :: tr_desc(NTR) !< Descriptions and metadata for the tracers in this package
end type ISOMIP_tracer_CS

contains


!> This subroutine is used to register tracer fields
function register_ISOMIP_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),      intent(in) :: HI    !<A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the open file
                                                       !! to parse for model parameter values.
  type(ISOMIP_tracer_CS),     pointer    :: CS !<A pointer that is set to point to the control
                                                       !! structure for this module (in/out).
  type(tracer_registry_type), pointer    :: tr_Reg !<A pointer to the tracer registry.
  type(MOM_restart_CS),       pointer    :: restart_CS !<A pointer to the restart control structure.

  character(len=80)  :: name, longname
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "ISOMIP_tracer" ! This module's name.
  character(len=200) :: inputdir
  character(len=48)  :: flux_units ! The units for tracer fluxes, usually
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_ISOMIP_tracer
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "ISOMIP_register_tracer called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "ISOMIP_TRACER_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial "//&
                 "conditions for the ISOMIP tracers, or blank to initialize "//&
                 "them internally.", default=" ")
  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    CS%tracer_IC_file = trim(inputdir)//trim(CS%tracer_IC_file)
    call log_param(param_file, mdl, "INPUTDIR/ISOMIP_TRACER_IC_FILE", &
                   CS%tracer_IC_file)
  endif
  call get_param(param_file, mdl, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified from MOM_initialization.F90.", default=.false.)

  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR)) ; CS%tr(:,:,:,:) = 0.0

  do m=1,NTR
    if (m < 10) then ; write(name,'("tr_D",I1.1)') m
    else ; write(name,'("tr_D",I2.2)') m ; endif
    write(longname,'("Concentration of ISOMIP Tracer ",I2.2)') m
    CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)
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
          flux_type=' ', implementation=' ', caller="register_ISOMIP_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  register_ISOMIP_tracer = .true.
end function register_ISOMIP_tracer

!> Initializes the NTR tracer fields in tr(:,:,:,:)
! and it sets up the tracer output.
subroutine initialize_ISOMIP_tracer(restart, day, G, GV, h, diag, OBC, CS, &
                                    ALE_sponge_CSp)

  type(ocean_grid_type),                 intent(in) :: G !< Grid structure.
  type(verticalGrid_type),               intent(in) :: GV !< The ocean's vertical grid structure.
  logical,                               intent(in) :: restart !< .true. if the fields have already
                                                       !! been read from a restart file.
  type(time_type), target,               intent(in) :: day !< Time of the start of the run.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< Layer thickness [H ~> m or kg m-2].
  type(diag_ctrl),               target, intent(in) :: diag !< A structure that is used to regulate
                                                       !! diagnostic output.
  type(ocean_OBC_type),                  pointer    :: OBC !< This open boundary condition type specifies
                                                       !! whether, where, and what open boundary conditions
                                                       !! are used. This is not being used for now.
  type(ISOMIP_tracer_CS),                pointer    :: CS !< The control structure returned by a previous call
                                                       !! to ISOMIP_register_tracer.
  type(ALE_sponge_CS),                   pointer    :: ALE_sponge_CSp !< A pointer to the control structure for
                                                       !! the sponges, if they are in use.  Otherwise this
                                                       !! may be unassociated.

  real, allocatable :: temp(:,:,:)
  real, pointer, dimension(:,:,:) :: &
    OBC_tr1_u => NULL(), & ! These arrays should be allocated and set to
    OBC_tr1_v => NULL()    ! specify the values of tracer 1 that should come
                           ! in through u- and v- points through the open
                           ! boundary conditions, in the same units as tr.
  character(len=16) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for tracer fluxes, usually
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  real, pointer :: tr_ptr(:,:,:) => NULL()
  real :: h_neglect         ! A thickness that is so small it is usually lost
                            ! in roundoff and can be neglected [H ~> m or kg m-2].
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB

  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  h_neglect = GV%H_subroundoff

  CS%Time => day
  CS%diag => diag

  if (.not.restart) then
    if (len_trim(CS%tracer_IC_file) >= 1) then
      !  Read the tracer concentrations from a netcdf file.
      if (.not.file_exists(CS%tracer_IC_file, G%Domain)) &
        call MOM_error(FATAL, "ISOMIP_initialize_tracer: Unable to open "// &
                        CS%tracer_IC_file)
      do m=1,NTR
        call query_vardesc(CS%tr_desc(m), name, caller="initialize_ISOMIP_tracer")
        call MOM_read_data(CS%tracer_IC_file, trim(name), CS%tr(:,:,:,m), G%Domain)
      enddo
    else
      do m=1,NTR
        do k=1,nz ; do j=js,je ; do i=is,ie
          CS%tr(i,j,k,m) = 0.0
        enddo ; enddo ; enddo
      enddo
    endif
  endif ! restart

! the following does not work in layer mode yet
!!  if ( CS%use_sponge ) then
  !   If sponges are used, this example damps tracers in sponges in the
  ! northern half of the domain to 1 and tracers in the southern half
  ! to 0.  For any tracers that are not damped in the sponge, the call
  ! to set_up_sponge_field can simply be omitted.
!    if (.not.associated(ALE_sponge_CSp)) &
!      call MOM_error(FATAL, "ISOMIP_initialize_tracer: "// &
!        "The pointer to ALEsponge_CSp must be associated if SPONGE is defined.")

!    allocate(temp(G%isd:G%ied,G%jsd:G%jed,nz))

!    do j=js,je ; do i=is,ie
!      if (G%geoLonT(i,j) >= 790.0 .AND. G%geoLonT(i,j) <= 800.0) then
!        temp(i,j,:) = 1.0
!      else
!        temp(i,j,:) = 0.0
!      endif
!    enddo ; enddo

      !   do m=1,NTR
!    do m=1,1
      ! This is needed to force the compiler not to do a copy in the sponge
      ! calls.  Curses on the designers and implementers of Fortran90.
!      tr_ptr => CS%tr(:,:,:,m)
!      call set_up_ALE_sponge_field(temp, G, tr_ptr, ALE_sponge_CSp)
!    enddo
!    deallocate(temp)
!  endif

end subroutine initialize_ISOMIP_tracer

!> This subroutine applies diapycnal diffusion, including the surface boundary
!! conditions and any other column tracer physics or chemistry to the tracers from this file.
subroutine ISOMIP_tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, GV, US, CS, &
                                        evap_CFL_limit, minimum_forcing_depth)
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
  type(ISOMIP_tracer_CS),  pointer    :: CS !< The control structure returned by a previous
                                              !! call to ISOMIP_register_tracer.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]

! The arguments to this subroutine are redundant in that
!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  real :: melt(SZI_(G),SZJ_(G)) ! melt water (positive for melting, negative for freezing) [R Z T-1 ~> kg m-2 s-1]
  real :: mmax                ! The global maximum melting rate [R Z T-1 ~> kg m-2 s-1]
  character(len=256) :: mesg  ! The text of an error message
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return

  melt(:,:) = fluxes%iceshelf_melt(:,:)

  ! max. melt
  mmax = MAXVAL(melt(is:ie,js:je))
  call max_across_PEs(mmax)
  ! write(mesg,*) 'max melt = ', mmax
  ! call MOM_mesg(mesg, 5)
  ! dye melt water (m=1), dye = 1 if melt=max(melt)
  do m=1,NTR
    do j=js,je ; do i=is,ie
      if (melt(i,j) > 0.0) then ! melting
        CS%tr(i,j,1:2,m) = melt(i,j)/mmax ! inject dye in the ML
      else ! freezing
        CS%tr(i,j,1:2,m) = 0.0
      endif
    enddo ; enddo
  enddo

  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,NTR
      do k=1,nz ;do j=js,je ; do i=is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m), dt, fluxes, h_work, &
                                          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  else
    do m=1,NTR
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

end subroutine ISOMIP_tracer_column_physics

!> This subroutine extracts the surface fields from this tracer package that
!! are to be shared with the atmosphere in coupled configurations.
!! This particular tracer package does not report anything back to the coupler.
subroutine ISOMIP_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(ISOMIP_tracer_CS),  pointer       :: CS !< The control structure returned by a previous
                                               !! call to ISOMIP_register_tracer.

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

end subroutine ISOMIP_tracer_surface_state

!> Deallocate any memory used by the ISOMIP tracer package
subroutine ISOMIP_tracer_end(CS)
  type(ISOMIP_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                        !! call to ISOMIP_register_tracer.
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end subroutine ISOMIP_tracer_end

end module ISOMIP_tracer
