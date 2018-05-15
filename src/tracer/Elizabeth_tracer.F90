!> This module contains the routines used to set up and use a set of (one for now)
!! dynamically passive tracers. For now, three passive tracers can be injected in
!! the domain
!! Set up and use passive tracers requires the following:
!! (1) register_Elizabeth_tracer
!! (2) apply diffusion, physics/chemistry and advect the tracer

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, 2002                                           *
!*  Adapted to the IDEAL_IS test case by Gustavo Marques, Oct 2016     *
!*********+*********+*********+*********+*********+*********+***********

module Elizabeth_tracer

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl
use MOM_diag_to_Z, only : register_Z_tracer, diag_to_Z_CS
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_hor_index, only : hor_index_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_restart, only : register_restart_field, MOM_restart_CS
use MOM_ALE_sponge, only : set_up_ALE_sponge_field, ALE_sponge_CS
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type, get_time
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_registry, only : add_tracer_diagnostics, add_tracer_OBC_values
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_variables, only : surface
use MOM_open_boundary, only : ocean_OBC_type
use MOM_verticalGrid, only : verticalGrid_type
use MOM_coms, only : max_across_PEs, min_across_PEs

implicit none ; private

#include <MOM_memory.h>

!< Publicly available functions
public register_Elizabeth_tracer, initialize_Elizabeth_tracer
public Elizabeth_tracer_column_physics, Elizabeth_tracer_end

!< ntr is the number of tracers in this module. (originally 2)
integer, parameter :: NTR = 1

type p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d

!> tracer control structure
type, public :: Elizabeth_tracer_CS ; private
  logical :: coupled_tracers = .false.  !< These tracers are not offered to the
                                        !< coupler.
  character(len = 200) :: tracer_IC_file !< The full path to the IC file, or " "
                                   !< to initialize internally.
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL()
  real, pointer :: tr(:,:,:,:) => NULL()   !< The array of tracers used in this
                                           !< subroutine, in g m-3?
  real, pointer :: tr_aux(:,:,:,:) => NULL() !< The masked tracer concentration
                                             !< for output, in g m-3.
  type(p3d), dimension(NTR) :: &
    tr_adx, &!< Tracer zonal advective fluxes in g m-3 m3 s-1.
    tr_ady, &!< Tracer meridional advective fluxes in g m-3 m3 s-1.
    tr_dfx, &!< Tracer zonal diffusive fluxes in g m-3 m3 s-1.
    tr_dfy   !< Tracer meridional diffusive fluxes in g m-3 m3 s-1.
  real :: land_val(NTR) = -1.0 !< The value of tr used where land is masked out.
  real :: lenlat           ! the latitudinal or y-direction length of the domain.
  real :: lenlon           ! the longitudinal or x-direction length of the domain.
  real :: CSL              ! The length of the continental shelf (x dir, km)
  real :: lensponge        ! the length of the sponge layer.
  logical :: mask_tracers  !< If true, tracers are masked out in massless layers.
  logical :: use_sponge

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                             !< timing of diagnostic output.
  integer, dimension(NTR) :: id_tracer = -1, id_tr_adx = -1, id_tr_ady = -1
  integer, dimension(NTR) :: id_tr_dfx = -1, id_tr_dfy = -1

  type(vardesc) :: tr_desc(NTR)
end type Elizabeth_tracer_CS

contains


!> This subroutine is used to register tracer fields
function register_Elizabeth_tracer(HI, GV, param_file, CS, tr_Reg, &
                                      restart_CS)
  type(hor_index_type),      intent(in) :: HI    !<A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !<A structure indicating the open file to parse for model parameter values.
  type(Elizabeth_tracer_CS),       pointer    :: CS !<A pointer that is set to point to the control structure for this module (in/out).
  type(tracer_registry_type), pointer    :: tr_Reg !<A pointer to the tracer registry.
  type(MOM_restart_CS),       pointer    :: restart_CS !<A pointer to the restart control structure.

  character(len=80)  :: name, longname
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "Elizabeth_tracer" ! This module's name.
  character(len=200) :: inputdir
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_Elizabeth_tracer
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "Elizabeth_register_tracer called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "Elizabeth_TRACER_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial \n"//&
                 "conditions for the Elizabeth tracers, or blank to initialize \n"//&
                 "them internally.", default=" ")
  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    CS%tracer_IC_file = trim(inputdir)//trim(CS%tracer_IC_file)
    call log_param(param_file, mod, "INPUTDIR/Elizabeth_TRACER_IC_FILE", &
                   CS%tracer_IC_file)
  endif
  call get_param(param_file, mod, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. \n"//&
                 "The exact location and properties of those sponges are \n"//&
                 "specified from MOM_initialization.F90.", default=.false.)

  call get_param(param_file, mod, "LENLAT", CS%lenlat, &
                 "The latitudinal or y-direction length of the domain", &
                 fail_if_missing=.true., do_not_log=.true.)

  call get_param(param_file, mod, "LENLON", CS%lenlon, &
                 "The longitudinal or x-direction length of the domain", &
                 fail_if_missing=.true., do_not_log=.true.)

  call get_param(param_file, mod, "CONT_SHELF_LENGTH", CS%CSL, &
                 "The length of the continental shelf (x dir, km).", &
                 default=15.0)

  call get_param(param_file, mod, "LENSPONGE", CS%lensponge, &
                 "The length of the sponge layer (km).", &
                 default=10.0)

  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR)) ; CS%tr(:,:,:,:) = 0.0
  if (CS%mask_tracers) then
    allocate(CS%tr_aux(isd:ied,jsd:jed,nz,NTR)) ; CS%tr_aux(:,:,:,:) = 0.0
  endif

  do m=1,NTR
    if (m < 10) then ; write(name,'("tr_D",I1.1)') m
    else ; write(name,'("tr_D",I2.2)') m ; endif
    write(longname,'("Concentration of Elizabeth Tracer ",I2.2)') m
    CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mod)

    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for the restart file.
    call register_restart_field(tr_ptr, CS%tr_desc(m), .true., restart_CS)
    ! Register the tracer for horizontal advection & diffusion.
    call register_tracer(tr_ptr, CS%tr_desc(m), param_file, HI, GV, tr_Reg, &
                         tr_desc_ptr=CS%tr_desc(m))

  enddo

  CS%tr_Reg => tr_Reg
  register_Elizabeth_tracer = .true.
end function register_Elizabeth_tracer

!> Initializes the NTR tracer fields in tr(:,:,:,:)
! and it sets up the tracer output.
subroutine initialize_Elizabeth_tracer(restart, day, G, GV, h, diag, OBC, CS, &
                                    sponge_CSp, diag_to_Z_CSp)

  type(ocean_grid_type),                 intent(in) :: G !< Grid structure.
  type(verticalGrid_type),               intent(in) :: GV !< The ocean's vertical grid structure.
  logical,                               intent(in) :: restart !< .true. if the fields have already been read from a restart file.
  type(time_type), target,               intent(in) :: day !< Time of the start of the run.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h !< Layer thickness, in m or kg m-2.
  type(diag_ctrl), target,               intent(in) :: diag
  type(ocean_OBC_type),                  pointer    :: OBC !< This open boundary condition type specifies whether, where, and what open boundary conditions are used. This is not being used for now.
  type(Elizabeth_tracer_CS),                pointer    :: CS !< The control structure returned by a previous call to Elizabeth_register_tracer.
  type(sponge_CS),                   pointer    :: sponge_CSp !< A pointer to the control structure for the sponges, if they are in use.  Otherwise this may be unassociated.
  type(diag_to_Z_CS),                    pointer    :: diag_to_Z_CSp !< A pointer to the control structure for diagnostics in depth space.

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
  real :: dist2  ! The distance squared from a line, in m2.
  real :: h_neglect         ! A thickness that is so small it is usually lost
                            ! in roundoff and can be neglected, in m.
  real :: e(SZK_(G)+1), e_top, e_bot, d_tr
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB

  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  h_neglect = GV%H_subroundoff

  CS%Time => day
  CS%diag => diag

  if (.not.restart) then
    if (len_trim(CS%tracer_IC_file) >= 1) then
      !  Read the tracer concentrations from a netcdf file.
      if (.not.file_exists(CS%tracer_IC_file, G%Domain)) &
        call MOM_error(FATAL, "Elizabeth_initialize_tracer: Unable to open "// &
                        CS%tracer_IC_file)
      do m=1,NTR
        call query_vardesc(CS%tr_desc(m), name, caller="initialize_Elizabeth_tracer")
        call read_data(CS%tracer_IC_file, trim(name), &
                       CS%tr(:,:,:,m), domain=G%Domain%mpp_domain)
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
                       caller="initialize_Elizabeth_tracer")
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

end subroutine initialize_Elizabeth_tracer

!> This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! This is a simple example of a set of advected passive tracers.
subroutine Elizabeth_tracer_column_physics(h_old, h_new,  ea,  eb, dt, G, GV, CS)
  type(ocean_grid_type),                 intent(in) :: G
  type(verticalGrid_type),               intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h_old, h_new, ea, eb
  real,                                  intent(in) :: dt
  type(Elizabeth_tracer_CS),                  pointer    :: CS

! Arguments: h_old -  Layer thickness before entrainment, in m or kg m-2.
!  (in)      h_new -  Layer thickness after entrainment, in m or kg m-2.
!  (in)      ea - an array to which the amount of fluid entrained
!                 from the layer above during this call will be
!                 added, in m or kg m-2.
!  (in)      eb - an array to which the amount of fluid entrained
!                 from the layer below during this call will be
!                 added, in m or kg m-2.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 Elizabeth_register_tracer.
!
! The arguments to this subroutine are redundant in that
!     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]

  real :: mmax, area
  real :: b1(SZI_(G))          ! b1 and c1 are variables used by the
  real :: c1(SZI_(G),SZK_(G))  ! tridiagonal solver.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_work ! Used so that h can be modified
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(CS)) return

  m=1
  do j=js,je ; do i=is,ie
    ! Set tracer to 1.0 in the surface of the continental shelf
    if (G%geoLonT(i,j) <= (CS%CSL) then
       CS%tr(i,j,1,m) = 1.0 ! first layer
    endif

  enddo ; enddo

  do j=js,je ; do i=is,ie
     ! remove tracer in the sponge layer
     if (G%geoLonT(i,j) >= (CS%lenlon - CS%lensponge) .AND. G%geoLonT(i,j) <= CS%lenlon) then
        CS%tr(i,j,:,m) = 0.0 ! all layers
     endif

  enddo ; enddo

  if (CS%mask_tracers) then
    do m = 1,NTR ; if (CS%id_tracer(m) > 0) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        if (h_new(i,j,k) < 1.1*GV%Angstrom) then
          CS%tr_aux(i,j,k,m) = CS%land_val(m)
        else
          CS%tr_aux(i,j,k,m) = CS%tr(i,j,k,m)
        endif
      enddo ; enddo ; enddo
    endif ; enddo
  endif

  do m=1,NTR
    if (CS%mask_tracers) then
      if (CS%id_tracer(m)>0) &
        call post_data(CS%id_tracer(m),CS%tr_aux(:,:,:,m),CS%diag)
    else
      if (CS%id_tracer(m)>0) &
        call post_data(CS%id_tracer(m),CS%tr(:,:,:,m),CS%diag)
    endif
    if (CS%id_tr_adx(m)>0) &
      call post_data(CS%id_tr_adx(m),CS%tr_adx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_ady(m)>0) &
      call post_data(CS%id_tr_ady(m),CS%tr_ady(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfx(m)>0) &
      call post_data(CS%id_tr_dfx(m),CS%tr_dfx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfy(m)>0) &
      call post_data(CS%id_tr_dfy(m),CS%tr_dfy(m)%p(:,:,:),CS%diag)
  enddo

end subroutine Elizabeth_tracer_column_physics

subroutine Elizabeth_tracer_end(CS)
  type(Elizabeth_tracer_CS), pointer :: CS
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    if (associated(CS%tr_aux)) deallocate(CS%tr_aux)
    do m=1,NTR
      if (associated(CS%tr_adx(m)%p)) deallocate(CS%tr_adx(m)%p)
      if (associated(CS%tr_ady(m)%p)) deallocate(CS%tr_ady(m)%p)
      if (associated(CS%tr_dfx(m)%p)) deallocate(CS%tr_dfx(m)%p)
      if (associated(CS%tr_dfy(m)%p)) deallocate(CS%tr_dfy(m)%p)
    enddo

    deallocate(CS)
  endif
end subroutine Elizabeth_tracer_end

end module Elizabeth_tracer
