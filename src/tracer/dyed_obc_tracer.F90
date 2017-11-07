module dyed_obc_tracer

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator,      only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,      only : diag_ctrl
use MOM_diag_to_Z,          only : register_Z_tracer, diag_to_Z_CS
use MOM_error_handler,      only : MOM_error, FATAL, WARNING
use MOM_file_parser,        only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,       only : forcing
use MOM_hor_index,          only : hor_index_type
use MOM_grid,               only : ocean_grid_type
use MOM_io,                 only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary,      only : ocean_OBC_type
use MOM_restart,            only : register_restart_field, MOM_restart_CS
use MOM_time_manager,       only : time_type, get_time
use MOM_tracer_registry,    only : register_tracer, tracer_registry_type
use MOM_tracer_registry,    only : add_tracer_diagnostics, add_tracer_OBC_values
use MOM_tracer_diabatic,    only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_variables,          only : surface
use MOM_verticalGrid,       only : verticalGrid_type

use coupler_types_mod,      only : coupler_type_set_data, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

public register_dyed_obc_tracer, initialize_dyed_obc_tracer
public dyed_obc_tracer_column_physics, dyed_obc_tracer_end

! ntr is the number of tracers in this module.
integer, parameter :: NTR = 4

type p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d

type, public :: dyed_obc_tracer_CS ; private
  logical :: coupled_tracers = .false.  ! These tracers are not offered to the
                                        ! coupler.
  character(len=200) :: tracer_IC_file ! The full path to the IC file, or " "
                                   ! to initialize internally.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL()
  real, pointer :: tr(:,:,:,:) => NULL()   ! The array of tracers used in this
                                           ! subroutine, in g m-3?
  type(p3d), dimension(NTR) :: &
    tr_adx, &! Tracer zonal advective fluxes in g m-3 m3 s-1.
    tr_ady, &! Tracer meridional advective fluxes in g m-3 m3 s-1.
    tr_dfx, &! Tracer zonal diffusive fluxes in g m-3 m3 s-1.
    tr_dfy   ! Tracer meridional diffusive fluxes in g m-3 m3 s-1.
  real :: land_val(NTR) = -1.0 ! The value of tr used where land is masked out.

  integer, dimension(NTR) :: ind_tr ! Indices returned by aof_set_coupler_flux
             ! if it is used and the surface tracer concentrations are to be
             ! provided to the coupler.

  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  integer, dimension(NTR) :: id_tracer = -1, id_tr_adx = -1, id_tr_ady = -1
  integer, dimension(NTR) :: id_tr_dfx = -1, id_tr_dfy = -1

  type(vardesc) :: tr_desc(NTR)
end type dyed_obc_tracer_CS

contains

!> This subroutine is used to register tracer fields and subroutines
!! to be used with MOM.
function register_dyed_obc_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(dyed_obc_tracer_CS),   pointer    :: CS   !< A pointer that is set to point to the
                                                 !! control structure for this module
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer to the tracer registry.
  type(MOM_restart_CS),       pointer    :: restart_CS !< A pointer to the restart control structure.

! Local variables
  character(len=80)  :: name, longname
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "dyed_obc_tracer" ! This module's name.
  character(len=200) :: inputdir
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_dyed_obc_tracer
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "dyed_obc_register_tracer called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "dyed_obc_TRACER_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial \n"//&
                 "conditions for the dyed_obc tracers, or blank to initialize \n"//&
                 "them internally.", default=" ")
  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    CS%tracer_IC_file = trim(inputdir)//trim(CS%tracer_IC_file)
    call log_param(param_file, mdl, "INPUTDIR/dyed_obc_TRACER_IC_FILE", &
                   CS%tracer_IC_file)
  endif

  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR)) ; CS%tr(:,:,:,:) = 0.0

  do m=1,NTR
    write(name,'("dye_",I1.1)') m
    write(longname,'("Concentration of dyed_obc Tracer ",I1.1)') m
    CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)

    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for the restart file.
    call register_restart_field(tr_ptr, CS%tr_desc(m), .true., restart_CS)
    ! Register the tracer for horizontal advection & diffusion.
    call register_tracer(tr_ptr, CS%tr_desc(m), param_file, HI, GV, tr_Reg, &
                         tr_desc_ptr=CS%tr_desc(m))

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = aof_set_coupler_flux(trim(name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_dyed_obc_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  register_dyed_obc_tracer = .true.
end function register_dyed_obc_tracer

!> This subroutine initializes the NTR tracer fields in tr(:,:,:,:)
!! and it sets up the tracer output.
subroutine initialize_dyed_obc_tracer(restart, day, G, GV, h, diag, OBC, CS, &
                                  diag_to_Z_CSp)
  type(ocean_grid_type),                 intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),               intent(in) :: GV   !< The ocean's vertical grid structure
  logical,                               intent(in) :: restart !< .true. if the fields have already
                                                               !! been read from a restart file.
  type(time_type), target,               intent(in) :: day     !< Time of the start of the run.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h    !< Layer thicknesses, in H (usually m or kg m-2)
  type(diag_ctrl), target,               intent(in) :: diag    !< Structure used to regulate diagnostic output.
  type(ocean_OBC_type),                  pointer    :: OBC     !< Structure specifying open boundary options.
  type(dyed_obc_tracer_CS),              pointer    :: CS      !< The control structure returned by a previous
                                                               !! call to dyed_obc_register_tracer.
  type(diag_to_Z_CS),                    pointer    :: diag_to_Z_CSp !< A pointer to the control structure
                                                                     !! for diagnostics in depth space.

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
  real :: h_neglect         ! A thickness that is so small it is usually lost
                            ! in roundoff and can be neglected, in m.
  real :: e(SZK_(G)+1), e_top, e_bot, d_tr
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
        call MOM_error(FATAL, "dyed_obc_initialize_tracer: Unable to open "// &
                        CS%tracer_IC_file)
      do m=1,NTR
        call query_vardesc(CS%tr_desc(m), name, caller="initialize_dyed_obc_tracer")
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

  ! This needs to be changed if the units of tracer are changed above.
  if (GV%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
  else ; flux_units = "kg s-1" ; endif

  do m=1,NTR
    ! Register the tracer for the restart file.
    call query_vardesc(CS%tr_desc(m), name, units=units, longname=longname, &
                       caller="initialize_dyed_obc_tracer")
    CS%id_tracer(m) = register_diag_field("ocean_model", trim(name), CS%diag%axesTL, &
        day, trim(longname) , trim(units))
    CS%id_tr_adx(m) = register_diag_field("ocean_model", trim(name)//"_adx", &
        CS%diag%axesCuL, day, trim(longname)//" advective zonal flux" , &
        trim(flux_units))
    CS%id_tr_ady(m) = register_diag_field("ocean_model", trim(name)//"_ady", &
        CS%diag%axesCvL, day, trim(longname)//" advective meridional flux" , &
        trim(flux_units))
    CS%id_tr_dfx(m) = register_diag_field("ocean_model", trim(name)//"_dfx", &
        CS%diag%axesCuL, day, trim(longname)//" diffusive zonal flux" , &
        trim(flux_units))
    CS%id_tr_dfy(m) = register_diag_field("ocean_model", trim(name)//"_dfy", &
        CS%diag%axesCvL, day, trim(longname)//" diffusive zonal flux" , &
        trim(flux_units))
    if (CS%id_tr_adx(m) > 0) call safe_alloc_ptr(CS%tr_adx(m)%p,IsdB,IedB,jsd,jed,nz)
    if (CS%id_tr_ady(m) > 0) call safe_alloc_ptr(CS%tr_ady(m)%p,isd,ied,JsdB,JedB,nz)
    if (CS%id_tr_dfx(m) > 0) call safe_alloc_ptr(CS%tr_dfx(m)%p,IsdB,IedB,jsd,jed,nz)
    if (CS%id_tr_dfy(m) > 0) call safe_alloc_ptr(CS%tr_dfy(m)%p,isd,ied,JsdB,JedB,nz)

!    Register the tracer for horizontal advection & diffusion.
    if ((CS%id_tr_adx(m) > 0) .or. (CS%id_tr_ady(m) > 0) .or. &
        (CS%id_tr_dfx(m) > 0) .or. (CS%id_tr_dfy(m) > 0)) &
      call add_tracer_diagnostics(name, CS%tr_Reg, CS%tr_adx(m)%p, &
                                  CS%tr_ady(m)%p, CS%tr_dfx(m)%p, CS%tr_dfy(m)%p)

    call register_Z_tracer(CS%tr(:,:,:,m), trim(name), longname, units, &
                           day, G, diag_to_Z_CSp)
  enddo

end subroutine initialize_dyed_obc_tracer

!> This subroutine applies diapycnal diffusion and any other column
!! tracer physics or chemistry to the tracers from this file.
!! This is a simple example of a set of advected passive tracers.
!!
!! The arguments to this subroutine are redundant in that
!!     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]
subroutine dyed_obc_tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, GV, CS, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),                 intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),               intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h_old !< Layer thickness before entrainment,
                                                                !! in m or kg m-2.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h_new !< Layer thickness after entrainment,
                                                                !! in m or kg m-2.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: ea    !< an array to which the amount of
                                              !! fluid entrained from the layer above during this
                                              !! call will be added, in m or kg m-2.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: eb    !< an array to which the amount of
                                              !! fluid entrained from the layer below during this
                                              !! call will be added, in m or kg m-2.
  type(forcing),                         intent(in) :: fluxes !< A structure containing pointers to
                                              !! any possible forcing fields.  Unused fields have NULL ptrs.
  real,                                  intent(in) :: dt   !< The amount of time covered by this call, in s
  type(dyed_obc_tracer_CS),                  pointer    :: CS   !< The control structure returned by a previous
                                                            !! call to dyed_obc_register_tracer.
  real,                        optional,intent(in)  :: evap_CFL_limit
  real,                        optional,intent(in)  :: minimum_forcing_depth

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
      enddo ; enddo ; enddo;
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m) , dt, fluxes, h_work, &
          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  else
    do m=1,NTR
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

  do m=1,NTR
    if (CS%id_tracer(m)>0) &
      call post_data(CS%id_tracer(m),CS%tr(:,:,:,m),CS%diag)
    if (CS%id_tr_adx(m)>0) &
      call post_data(CS%id_tr_adx(m),CS%tr_adx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_ady(m)>0) &
      call post_data(CS%id_tr_ady(m),CS%tr_ady(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfx(m)>0) &
      call post_data(CS%id_tr_dfx(m),CS%tr_dfx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfy(m)>0) &
      call post_data(CS%id_tr_dfy(m),CS%tr_dfy(m)%p(:,:,:),CS%diag)
  enddo

end subroutine dyed_obc_tracer_column_physics

!> Clean up memory allocations, if any.
subroutine dyed_obc_tracer_end(CS)
  type(dyed_obc_tracer_CS), pointer :: CS
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    do m=1,NTR
      if (associated(CS%tr_adx(m)%p)) deallocate(CS%tr_adx(m)%p)
      if (associated(CS%tr_ady(m)%p)) deallocate(CS%tr_ady(m)%p)
      if (associated(CS%tr_dfx(m)%p)) deallocate(CS%tr_dfx(m)%p)
      if (associated(CS%tr_dfy(m)%p)) deallocate(CS%tr_dfy(m)%p)
    enddo

    deallocate(CS)
  endif
end subroutine dyed_obc_tracer_end

!> \namespace dyed_obc_tracer
!!                                                                     *
!!  By Kate Hedstrom, 2017, copied from DOME tracers.                  *
!!                                                                     *
!!    This file contains an example of the code that is needed to set  *
!!  up and use a set of dynamically passive tracers. These tracers     *
!!  dye the inflowing water, one per open boundary segment.            *
!!                                                                     *
!!    A single subroutine is called from within each file to register  *
!!  each of the tracers for reinitialization and advection and to      *
!!  register the subroutine that initializes the tracers and set up    *
!!  their output and the subroutine that does any tracer physics or    *
!!  chemistry along with diapycnal mixing (included here because some  *
!!  tracers may float or swim vertically or dye diapycnal processes).  *
!!                                                                     *
!!                                                                     *
!!  Macros written all in capital letters are defined in MOM_memory.h. *
!!                                                                     *
!!     A small fragment of the grid is shown below:                    *
!!                                                                     *
!!    j+1  x ^ x ^ x   At x:  q                                        *
!!    j+1  > o > o >   At ^:  v                                        *
!!    j    x ^ x ^ x   At >:  u                                        *
!!    j    > o > o >   At o:  h, tr                                    *
!!    j-1  x ^ x ^ x                                                   *
!!        i-1  i  i+1  At x & ^:                                       *
!!           i  i+1    At > & o:                                       *
!!                                                                     *
!!  The boundaries always run through q grid points (x).               *
!!                                                                     *
!!*******+*********+*********+*********+*********+*********+*********+**

end module dyed_obc_tracer
