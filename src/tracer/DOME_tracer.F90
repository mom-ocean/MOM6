!> A tracer package that is used as a diagnostic in the DOME experiments
module DOME_tracer

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator,   only : diag_ctrl
use MOM_error_handler,   only : MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_hor_index,       only : hor_index_type
use MOM_grid,            only : ocean_grid_type
use MOM_io,              only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary,   only : ocean_OBC_type, OBC_segment_tracer_type
use MOM_open_boundary,   only : OBC_segment_type
use MOM_restart,         only : MOM_restart_CS
use MOM_sponge,          only : set_up_sponge_field, sponge_CS
use MOM_time_manager,    only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface
use MOM_verticalGrid,    only : verticalGrid_type

use coupler_types_mod,      only : coupler_type_set_data, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

public register_DOME_tracer, initialize_DOME_tracer
public DOME_tracer_column_physics, DOME_tracer_surface_state, DOME_tracer_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

integer, parameter :: ntr = 11 !< The number of tracers in this module.

!> The DOME_tracer control structure
type, public :: DOME_tracer_CS ; private
  logical :: coupled_tracers = .false. !< These tracers are not offered to the coupler.
  character(len=200) :: tracer_IC_file !< The full path to the IC file, or " " to initialize internally.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL()   !< The array of tracers used in this package, in g m-3?
  real :: land_val(NTR) = -1.0 !< The value of tr used where land is masked out.
  logical :: use_sponge    !< If true, sponges may be applied somewhere in the domain.

  integer, dimension(NTR) :: ind_tr !< Indices returned by aof_set_coupler_flux if it is used and the
                                    !! surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.

  type(vardesc) :: tr_desc(NTR) !< Descriptions and metadata for the tracers
end type DOME_tracer_CS

contains

!> Register tracer fields and subroutines to be used with MOM.
function register_DOME_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(DOME_tracer_CS),       pointer    :: CS   !< A pointer that is set to point to the
                                                 !! control structure for this module
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer to the tracer registry.
  type(MOM_restart_CS),       pointer    :: restart_CS !< A pointer to the restart control structure.

! Local variables
  character(len=80)  :: name, longname
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "DOME_tracer" ! This module's name.
  character(len=48) :: flux_units ! The units for tracer fluxes, usually
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  character(len=200) :: inputdir
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_DOME_tracer
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "DOME_register_tracer called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DOME_TRACER_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial "//&
                 "conditions for the DOME tracers, or blank to initialize "//&
                 "them internally.", default=" ")
  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    CS%tracer_IC_file = trim(inputdir)//trim(CS%tracer_IC_file)
    call log_param(param_file, mdl, "INPUTDIR/DOME_TRACER_IC_FILE", &
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
    write(longname,'("Concentration of DOME Tracer ",I2.2)') m
    CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)
    if (GV%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
    else ; flux_units = "kg s-1" ; endif

    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                         name=name, longname=longname, units="kg kg-1", &
                         registry_diags=.true., restart_CS=restart_CS, &
                         flux_units=trim(flux_units), flux_scale=GV%H_to_MKS)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = aof_set_coupler_flux(trim(name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_DOME_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  register_DOME_tracer = .true.
end function register_DOME_tracer

!> Initializes the NTR tracer fields in tr(:,:,:,:) and sets up the tracer output.
subroutine initialize_DOME_tracer(restart, day, G, GV, US, h, diag, OBC, CS, &
                                  sponge_CSp, param_file)
  type(ocean_grid_type),                 intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),               intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                 intent(in) :: US   !< A dimensional unit scaling type
  logical,                               intent(in) :: restart !< .true. if the fields have already
                                                               !! been read from a restart file.
  type(time_type), target,               intent(in) :: day     !< Time of the start of the run.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl), target,               intent(in) :: diag    !< Structure used to regulate diagnostic output.
  type(ocean_OBC_type),                  pointer    :: OBC     !< Structure specifying open boundary options.
  type(DOME_tracer_CS),                  pointer    :: CS      !< The control structure returned by a previous
                                                               !! call to DOME_register_tracer.
  type(sponge_CS),                       pointer    :: sponge_CSp    !< A pointer to the control structure
                                                                     !! for the sponges, if they are in use.
  type(param_file_type),                 intent(in) :: param_file !< A structure to parse for run-time parameters

! Local variables
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
  real :: PI     ! 3.1415926... calculated as 4*atan(1)
  real :: tr_y   ! Initial zonally uniform tracer concentrations.
  real :: h_neglect         ! A thickness that is so small it is usually lost
                            ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: e(SZK_(G)+1), e_top, e_bot ! Heights [Z ~> m].
  real :: d_tr   ! A change in tracer concentraions, in tracer units.
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
        call MOM_error(FATAL, "DOME_initialize_tracer: Unable to open "// &
                        CS%tracer_IC_file)
      do m=1,NTR
        call query_vardesc(CS%tr_desc(m), name, caller="initialize_DOME_tracer")
        call MOM_read_data(CS%tracer_IC_file, trim(name), CS%tr(:,:,:,m), G%Domain)
      enddo
    else
      do m=1,NTR
        do k=1,nz ; do j=js,je ; do i=is,ie
          CS%tr(i,j,k,m) = 1.0e-20 ! This could just as well be 0.
        enddo ; enddo ; enddo
      enddo

!    This sets a stripe of tracer across the basin.
      do m=2,NTR ; do j=js,je ; do i=is,ie
        tr_y = 0.0
        if ((m <= 6) .and. (G%geoLatT(i,j) > (300.0+50.0*real(m-1))) .and. &
            (G%geoLatT(i,j) < (350.0+50.0*real(m-1)))) tr_y = 1.0
        do k=1,nz
!      This adds the stripes of tracer to every layer.
            CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + tr_y
        enddo
      enddo ; enddo ; enddo

      if (NTR > 7) then
        do j=js,je ; do i=is,ie
          e(nz+1) = -G%bathyT(i,j)
          do k=nz,1,-1
            e(K) = e(K+1) + h(i,j,k)*GV%H_to_Z
            do m=7,NTR
              e_top = (-600.0*real(m-1) + 3000.0) * US%m_to_Z
              e_bot = (-600.0*real(m-1) + 2700.0) * US%m_to_Z
              if (e_top < e(K)) then
                if (e_top < e(K+1)) then ; d_tr = 0.0
                elseif (e_bot < e(K+1)) then
                  d_tr = 1.0 * (e_top-e(K+1)) / ((h(i,j,k)+h_neglect)*GV%H_to_Z)
                else ; d_tr = 1.0 * (e_top-e_bot) / ((h(i,j,k)+h_neglect)*GV%H_to_Z)
                endif
              elseif (e_bot < e(K)) then
                if (e_bot < e(K+1)) then ; d_tr = 1.0
                else ; d_tr = 1.0 * (e(K)-e_bot) / ((h(i,j,k)+h_neglect)*GV%H_to_Z)
                endif
              else
                d_tr = 0.0
              endif
              if (h(i,j,k) < 2.0*GV%Angstrom_H) d_tr=0.0
              CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + d_tr
            enddo
          enddo
        enddo ; enddo
      endif

    endif
  endif ! restart

  if ( CS%use_sponge ) then
!   If sponges are used, this example damps tracers in sponges in the
! northern half of the domain to 1 and tracers in the southern half
! to 0.  For any tracers that are not damped in the sponge, the call
! to set_up_sponge_field can simply be omitted.
    if (.not.associated(sponge_CSp)) &
      call MOM_error(FATAL, "DOME_initialize_tracer: "// &
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
      call set_up_sponge_field(temp, tr_ptr, G, nz, sponge_CSp)
    enddo
    deallocate(temp)
  endif

end subroutine initialize_DOME_tracer

!> This subroutine applies diapycnal diffusion and any other column
!! tracer physics or chemistry to the tracers from this file.
!! This is a simple example of a set of advected passive tracers.
!!
!! The arguments to this subroutine are redundant in that
!!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)
subroutine DOME_tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, GV, US, CS, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(DOME_tracer_CS),    pointer    :: CS   !< The control structure returned by a previous
                                              !! call to DOME_register_tracer.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]

! Local variables
  real :: b1(SZI_(G))          ! b1 and c1 are variables used by the
  real :: c1(SZI_(G),SZK_(G))  ! tridiagonal solver.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_work ! Used so that h can be modified
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return

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

end subroutine DOME_tracer_column_physics

!> This subroutine extracts the surface fields from this tracer package that
!! are to be shared with the atmosphere in coupled configurations.
!! This particular tracer package does not report anything back to the coupler.
subroutine DOME_tracer_surface_state(state, h, G, CS)
  type(ocean_grid_type),  intent(in)    :: G  !< The ocean's grid structure.
  type(surface),          intent(inout) :: state !< A structure containing fields that
                                              !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                          intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(DOME_tracer_CS),   pointer       :: CS !< The control structure returned by a previous
                                              !! call to DOME_register_tracer.

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

  integer :: m, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,NTR
      !   This call loads the surface values into the appropriate array in the
      ! coupler-type structure.
      call coupler_type_set_data(CS%tr(:,:,1,m), CS%ind_tr(m), ind_csurf, &
                   state%tr_fields, idim=(/isd, is, ie, ied/), &
                   jdim=(/jsd, js, je, jed/) )
    enddo
  endif

end subroutine DOME_tracer_surface_state

!> Clean up memory allocations, if any.
subroutine DOME_tracer_end(CS)
  type(DOME_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                      !! call to DOME_register_tracer.
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end subroutine DOME_tracer_end

!> \namespace dome_tracer
!!
!!  By Robert Hallberg, 2002
!!
!!    This file contains an example of the code that is needed to set
!!  up and use a set (in this case eleven) of dynamically passive
!!  tracers.  These tracers dye the inflowing water or water initially
!!  within a range of latitudes or water initially in a range of
!!  depths.
!!
!!    A single subroutine is called from within each file to register
!!  each of the tracers for reinitialization and advection and to
!!  register the subroutine that initializes the tracers and set up
!!  their output and the subroutine that does any tracer physics or
!!  chemistry along with diapycnal mixing (included here because some
!!  tracers may float or swim vertically or dye diapycnal processes).

end module DOME_tracer
