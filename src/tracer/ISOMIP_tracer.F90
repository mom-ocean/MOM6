!> This module contains the routines used to set up and use a set of (one for now)
!! dynamically passive tracers. For now, just one passive tracer is injected in
!! the sponge layer.
!! Set up and use passive tracers requires the following:
!! (1) register_ISOMIP_tracer
!! (2)

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, 2002                                           *
!*  Adapted to the ISOMIP test case by Gustavo Marques, May 2016       *            !*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

module ISOMIP_tracer

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl
use MOM_diag_to_Z, only : register_Z_tracer, diag_to_Z_CS
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_hor_index, only : hor_index_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_restart, only : register_restart_field, MOM_restart_CS
use MOM_ALE_sponge, only : set_up_ALE_sponge_field, ALE_sponge_CS
use MOM_time_manager, only : time_type, get_time
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_registry, only : add_tracer_diagnostics, add_tracer_OBC_values
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_variables, only : surface
use MOM_open_boundary, only : ocean_OBC_type
use MOM_verticalGrid, only : verticalGrid_type
use MOM_coms, only : max_across_PEs
use coupler_util, only : set_coupler_values, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

!< Publicly available functions
public register_ISOMIP_tracer, initialize_ISOMIP_tracer
public ISOMIP_tracer_column_physics, ISOMIP_tracer_surface_state, ISOMIP_tracer_end

!< ntr is the number of tracers in this module.
integer, parameter :: ntr = 1

type p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d

!> tracer control structure
type, public :: ISOMIP_tracer_CS ; private
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
  logical :: mask_tracers  !< If true, tracers are masked out in massless layers.
  logical :: use_sponge

  integer, dimension(NTR) :: ind_tr !< Indices returned by aof_set_coupler_flux
             !< if it is used and the surface tracer concentrations are to be
             !< provided to the coupler.

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                             !< timing of diagnostic output.
  integer, dimension(NTR) :: id_tracer = -1, id_tr_adx = -1, id_tr_ady = -1
  integer, dimension(NTR) :: id_tr_dfx = -1, id_tr_dfy = -1

  type(vardesc) :: tr_desc(NTR)
end type ISOMIP_tracer_CS

contains


!> This subroutine is used to register tracer fields
function register_ISOMIP_tracer(HI, GV, param_file, CS, tr_Reg, &
                                      restart_CS)
  type(hor_index_type),      intent(in) :: HI    !<A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !<A structure indicating the open file to parse for model parameter values.
  type(ISOMIP_tracer_CS),       pointer    :: CS !<A pointer that is set to point to the control structure for this module (in/out).
  type(tracer_registry_type), pointer    :: tr_Reg !<A pointer to the tracer registry.
  type(MOM_restart_CS),       pointer    :: restart_CS !<A pointer to the restart control structure.

  character(len=80)  :: name, longname
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "ISOMIP_tracer" ! This module's name.
  character(len=200) :: inputdir
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
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "ISOMIP_TRACER_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial \n"//&
                 "conditions for the ISOMIP tracers, or blank to initialize \n"//&
                 "them internally.", default=" ")
  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    CS%tracer_IC_file = trim(inputdir)//trim(CS%tracer_IC_file)
    call log_param(param_file, mod, "INPUTDIR/ISOMIP_TRACER_IC_FILE", &
                   CS%tracer_IC_file)
  endif
  call get_param(param_file, mod, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. \n"//&
                 "The exact location and properties of those sponges are \n"//&
                 "specified from MOM_initialization.F90.", default=.false.)

  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR)) ; CS%tr(:,:,:,:) = 0.0
  if (CS%mask_tracers) then
    allocate(CS%tr_aux(isd:ied,jsd:jed,nz,NTR)) ; CS%tr_aux(:,:,:,:) = 0.0
  endif

  do m=1,NTR
    if (m < 10) then ; write(name,'("tr_D",I1.1)') m
    else ; write(name,'("tr_D",I2.2)') m ; endif
    write(longname,'("Concentration of ISOMIP Tracer ",I2.2)') m
    CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mod)

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
          flux_type=' ', implementation=' ', caller="register_ISOMIP_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  register_ISOMIP_tracer = .true.
end function register_ISOMIP_tracer

!> Initializes the NTR tracer fields in tr(:,:,:,:)
! and it sets up the tracer output.
subroutine initialize_ISOMIP_tracer(restart, day, G, GV, h, diag, OBC, CS, &
                                    ALE_sponge_CSp, diag_to_Z_CSp)

  type(ocean_grid_type),                 intent(in) :: G !< Grid structure.
  type(verticalGrid_type),               intent(in) :: GV !< The ocean's vertical grid structure.
  logical,                               intent(in) :: restart !< .true. if the fields have already been read from a restart file.
  type(time_type), target,               intent(in) :: day !< Time of the start of the run.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h !< Layer thickness, in m or kg m-2.
  type(diag_ctrl), target,               intent(in) :: diag
  type(ocean_OBC_type),                  pointer    :: OBC !< This open boundary condition type specifies whether, where, and what open boundary conditions are used. This is not being used for now.
  type(ISOMIP_tracer_CS),                pointer    :: CS !< The control structure returned by a previous call to ISOMIP_register_tracer.
  type(ALE_sponge_CS),                   pointer    :: ALE_sponge_CSp !< A pointer to the control structure for the sponges, if they are in use.  Otherwise this may be unassociated.
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
        call MOM_error(FATAL, "ISOMIP_initialize_tracer: Unable to open "// &
                        CS%tracer_IC_file)
      do m=1,NTR
        call query_vardesc(CS%tr_desc(m), name, caller="initialize_ISOMIP_tracer")
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

  ! This needs to be changed if the units of tracer are changed above.
  if (GV%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
  else ; flux_units = "kg s-1" ; endif

  do m=1,NTR
    ! Register the tracer for the restart file.
    call query_vardesc(CS%tr_desc(m), name, units=units, longname=longname, &
                       caller="initialize_ISOMIP_tracer")
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

end subroutine initialize_ISOMIP_tracer

!> This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! This is a simple example of a set of advected passive tracers.
subroutine ISOMIP_tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, GV, CS, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),                 intent(in) :: G
  type(verticalGrid_type),               intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h_old, h_new, ea, eb
  type(forcing),                         intent(in) :: fluxes
  real,                                  intent(in) :: dt
  type(ISOMIP_tracer_CS),                  pointer    :: CS
  real,                             optional,intent(in)  :: evap_CFL_limit
  real,                             optional,intent(in)  :: minimum_forcing_depth

! Arguments: h_old -  Layer thickness before entrainment, in m or kg m-2.
!  (in)      h_new -  Layer thickness after entrainment, in m or kg m-2.
!  (in)      ea - an array to which the amount of fluid entrained
!                 from the layer above during this call will be
!                 added, in m or kg m-2.
!  (in)      eb - an array to which the amount of fluid entrained
!                 from the layer below during this call will be
!                 added, in m or kg m-2.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 ISOMIP_register_tracer.
!
! The arguments to this subroutine are redundant in that
!     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]

  real :: mmax
  real :: b1(SZI_(G))          ! b1 and c1 are variables used by the
  real :: c1(SZI_(G),SZK_(G))  ! tridiagonal solver.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_work ! Used so that h can be modified
  real :: melt(SZI_(G),SZJ_(G))  ! melt water (positive for melting
                                 ! negative for freezing)
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(CS)) return

  melt(:,:) = fluxes%iceshelf_melt

  ! max. melt
  mmax = MAXVAL(melt(is:ie,js:je))
  call max_across_PEs(mmax)
  !write(*,*)'max melt', mmax
  ! dye melt water (m=1), dye = 1 if melt=max(melt)
  do m=1,NTR
     do j=js,je ; do i=is,ie
      if (melt(i,j) > 0.0) then ! melting
         !write(*,*)'i,j,melt,melt/mmax',i,j,melt(i,j),melt(i,j)/mmax
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

end subroutine ISOMIP_tracer_column_physics

!> This particular tracer package does not report anything back to the coupler.
! The code that is here is just a rough guide for packages that would.
subroutine ISOMIP_tracer_surface_state(state, h, G, CS)
  type(ocean_grid_type),                    intent(in)    :: G !< The ocean's grid structure.
  type(surface),                            intent(inout) :: state !< A structure containing fields that describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h !< Layer thickness, in m or kg m-2.
  type(ISOMIP_tracer_CS),                     pointer       :: CS !< The control structure returned by a previous call to ISOMIP_register_tracer.
  integer :: i, j, m, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,ntr
      !   This call loads the surface vlues into the appropriate array in the
      ! coupler-type structure.
      call set_coupler_values(CS%tr(:,:,1,1), state%tr_fields, CS%ind_tr(m), &
                              ind_csurf, is, ie, js, je)
    enddo
  endif

end subroutine ISOMIP_tracer_surface_state

subroutine ISOMIP_tracer_end(CS)
  type(ISOMIP_tracer_CS), pointer :: CS
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
end subroutine ISOMIP_tracer_end

end module ISOMIP_tracer
