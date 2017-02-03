!> Implements a boundary impulse response tracer to calculate Green's functions
module boundary_impulse_tracer
! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl
use MOM_diag_to_Z, only : register_Z_tracer, diag_to_Z_CS
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_io, only : file_exists, read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : register_restart_field, query_initialized, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type, get_time
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_registry, only : add_tracer_diagnostics, add_tracer_OBC_values
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init, only : tracer_Z_init
use MOM_variables, only : surface
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use coupler_util, only : set_coupler_values, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

public register_boundary_impulse_tracer, initialize_boundary_impulse_tracer
public boundary_impulse_tracer_column_physics, boundary_impulse_tracer_surface_state
public boundary_impulse_stock, boundary_impulse_tracer_end

! NTR_MAX is the maximum number of tracers in this module.
integer, parameter :: NTR_MAX = 1

type p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d

type, public :: boundary_impulse_tracer_CS ; private
  integer :: ntr=NTR_MAX    ! The number of tracers that are actually used.
  logical :: coupled_tracers = .false.  ! These tracers are not offered to the
                                        ! coupler.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL()
  real, pointer :: tr(:,:,:,:) => NULL()   ! The array of tracers used in this
                                           ! subroutine, in g m-3?
  type(p3d), dimension(NTR_MAX) :: &
    tr_adx, &! Tracer zonal advective fluxes in g m-3 m3 s-1.An Error Has Occurred


    tr_ady, &! Tracer meridional advective fluxes in g m-3 m3 s-1.
    tr_dfx, &! Tracer zonal diffusive fluxes in g m-3 m3 s-1.
    tr_dfy   ! Tracer meridional diffusive fluxes in g m-3 m3 s-1.
  logical :: mask_tracers  ! If true, boundary_impulse is masked out in massless layers.
  logical :: tracers_may_reinit  ! If true, boundary_impulse can be initialized if
                                 ! not found in restart file
  integer, dimension(NTR_MAX) :: &
    ind_tr, &  ! Indices returned by aof_set_coupler_flux if it is used and the
               ! surface tracer concentrations are to be provided to the coupler.
    id_tracer = -1, id_tr_adx = -1, id_tr_ady = -1, &
    id_tr_dfx = -1, id_tr_dfy = -1

  integer :: nkml ! Number of layers in mixed layer
  real, dimension(NTR_MAX)  :: land_val = -1.0
  real :: kw_eff ! An effective piston velocity used to flux tracer out at the surface
  real :: remaining_source_time ! How much longer (same units as the timestep) to
                                ! inject the tracer at the surface

  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()

  type(vardesc) :: tr_desc(NTR_MAX)
end type boundary_impulse_tracer_CS

contains

!> Read in runtime options and add boundary impulse tracer to tracer registry
function register_boundary_impulse_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI
  type(verticalGrid_type),    intent(in) :: GV
  type(param_file_type),      intent(in) :: param_file
  type(boundary_impulse_tracer_CS),        pointer    :: CS
  type(tracer_registry_type), pointer    :: tr_Reg
  type(MOM_restart_CS),       pointer    :: restart_CS
! This subroutine is used to register tracer fields and subroutines
! to be used with MOM.
! Arguments: HI - A horizontal index type structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in/out)  tr_Reg - A pointer that is set to point to the control structure
!                  for the tracer advection and diffusion module.
!  (in)      restart_CS - A pointer to the restart control structure.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "boundary_impulse_tracer" ! This module's name.
  character(len=200) :: inputdir ! The directory where the input files are.
  character(len=48)  :: var_name ! The variable's name.
  character(len=3)   :: name_tag ! String for creating identifying boundary_impulse
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_boundary_impulse_tracer
  integer :: isd, ied, jsd, jed, nz, m, i, j
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "register_boundary_impulse_tracer called with an "// &
                             "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "IMPULSE_SOURCE_TIME", CS%remaining_source_time, &
                 "Length of time for the boundary tracer to be injected\n"//&
                 "into the mixed layer. After this time has elapsed, the\n"//&
                 "surface becomes a sink for the boundary impulse tracer.", &
                 default=31536000.0)

  CS%ntr = NTR_MAX
  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr)) ; CS%tr(:,:,:,:) = 0.0

  CS%nkml = max(GV%nkml,1)

  do m=1,CS%ntr
    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    CS%tr_desc(m) = var_desc(trim("boundary_impulse"), "kg", &
        "Boundary impulse tracer", caller=mod)
    tr_ptr => CS%tr(:,:,:,m)
    call query_vardesc(CS%tr_desc(m), name=var_name, caller="register_boundary_impulse_tracer")
    ! Register the tracer for the restart file.
    call register_restart_field(tr_ptr, CS%tr_desc(m), &
                                .not. CS%tracers_may_reinit, restart_CS)
    ! Register the tracer for horizontal advection & diffusion.
    call register_tracer(tr_ptr, CS%tr_desc(m), param_file, HI, GV, tr_Reg, &
                         tr_desc_ptr=CS%tr_desc(m))

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = aof_set_coupler_flux(trim(var_name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_boundary_impulse_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_boundary_impulse_tracer = .true.

end function register_boundary_impulse_tracer

!> Initialize tracer from restart or set to 1 at surface to initialize
subroutine initialize_boundary_impulse_tracer(restart, day, G, GV, h, diag, OBC, CS, &
                                  sponge_CSp, diag_to_Z_CSp, tv)
  logical,                            intent(in) :: restart
  type(time_type), target,            intent(in) :: day
  type(ocean_grid_type),              intent(in) :: G
  type(verticalGrid_type),            intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
  type(diag_ctrl), target,            intent(in) :: diag
  type(ocean_OBC_type),               pointer    :: OBC
  type(boundary_impulse_tracer_CS),          pointer    :: CS
  type(sponge_CS),                    pointer    :: sponge_CSp
  type(diag_to_Z_CS),                 pointer    :: diag_to_Z_CSp
  type(thermo_var_ptrs),              intent(in) :: tv
!   This subroutine initializes the CS%ntr tracer fields in tr(:,:,:,:)
! and it sets up the tracer output.

! Arguments: restart - .true. if the fields have already been read from
!                     a restart file.
!  (in)      day - Time of the start of the run.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in)      OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in/out)  CS - The control structure returned by a previous call to
!                 register_boundary_impulse_tracer.
!  (in/out)  sponge_CSp - A pointer to the control structure for the sponges, if
!                         they are in use.  Otherwise this may be unassociated.
!  (in/out)  diag_to_Z_Csp - A pointer to the control structure for diagnostics
!                            in depth space.
  character(len=16) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for age tracer fluxes, either
                                ! years m3 s-1 or years kg s-1.
  logical :: OK
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%Time => day
  CS%diag => diag
  name = "boundary_impulse"

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=name, caller="initialize_boundary_impulse_tracer")
    if ((.not.restart) .or. (.not. &
        query_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp))) then
      do k=1,CS%nkml ; do j=jsd,jed ; do i=isd,ied
        CS%tr(i,j,k,m) = 1.0
      enddo ; enddo ; enddo
    endif
  enddo ! Tracer loop

  if (associated(OBC)) then
  ! All tracers but the first have 0 concentration in their inflows. As this
  ! is the default value, the following calls are unnecessary.
  ! do m=1,CS%ntr
  !  call add_tracer_OBC_values(trim(CS%tr_desc(m)%name), CS%tr_Reg, 0.0)
  ! enddo
  endif

  ! This needs to be changed if the units of tracer are changed above.
  if (GV%Boussinesq) then ; flux_units = "g salt/(m^2 s)"
  else ; flux_units = "g salt/(m^2 s)" ; endif

  do m=1,CS%ntr
    ! Register the tracer for the restart file.
    call query_vardesc(CS%tr_desc(m), name, units=units, longname=longname, &
                       caller="initialize_boundary_impulse_tracer")
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
                                  CS%tr_ady(m)%p,CS%tr_dfx(m)%p,CS%tr_dfy(m)%p)

    call register_Z_tracer(CS%tr(:,:,:,m), trim(name), longname, units, &
                           day, G, diag_to_Z_CSp)
  enddo

end subroutine initialize_boundary_impulse_tracer

! Apply source or sink at boundary and do vertical diffusion
subroutine boundary_impulse_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, CS, tv, debug, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),              intent(in) :: G
  type(verticalGrid_type),            intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h_old, h_new, ea, eb
  type(forcing),                      intent(in) :: fluxes
  real,                               intent(in) :: dt
  type(boundary_impulse_tracer_CS),                pointer    :: CS
  type(thermo_var_ptrs),              intent(in) :: tv
  logical,                            intent(in) :: debug
  real,                             optional,intent(in)  :: evap_CFL_limit
  real,                             optional,intent(in)  :: minimum_forcing_depth

!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! This is a simple example of a set of advected passive tracers.

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
!                 register_boundary_impulse_tracer.
!  (in)      tv - Thermodynamic structure with T and S
!  (in)      evap_CFL_limit - Limits how much water can be fluxed out of the top layer
!                             Stored previously in diabatic CS.
!  (in)      minimum_forcing_depth - The smallest depth over which fluxes can be applied
!                             Stored previously in diabatic CS.
!  (in)      debug - Calculates checksums
!
! The arguments to this subroutine are redundant in that
!     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]

  real :: Isecs_per_year = 1.0 / (365.0*86400.0)
  real :: year, h_total, scale, htot, Ih_limit
  integer :: secs, days
  integer :: i, j, k, is, ie, js, je, nz, m, k_max
  real, allocatable :: local_tr(:,:,:)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_work ! Used so that h can be modified

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  ! This uses applyTracerBoundaryFluxesInOut, usually in ALE mode
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do k=1,nz ;do j=js,je ; do i=is,ie
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo;
    call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,1), dt, fluxes, h_work, &
      evap_CFL_limit, minimum_forcing_depth)
    call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,1), G, GV)
  else
    call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,1), G, GV)
  endif

  ! Set surface conditions
  do m=1,1
    if(CS%remaining_source_time>0.0) then
      do k=1,CS%nkml ; do j=js,je ; do i=is,ie
        CS%tr(i,j,k,m) = 1.0
      enddo ; enddo ; enddo
      CS%remaining_source_time = CS%remaining_source_time-dt
    else
      do k=1,CS%nkml ; do j=js,je ; do i=is,ie
        CS%tr(i,j,k,m) = 0.0
      enddo ; enddo ; enddo
    endif

  enddo

  allocate(local_tr(G%isd:G%ied,G%jsd:G%jed,nz))
  do m=1,1
    if (CS%id_tracer(m)>0) then
      if (CS%mask_tracers) then
        do k=1,nz ; do j=js,je ; do i=is,ie
          if (h_new(i,j,k) < 1.1*GV%Angstrom) then
            local_tr(i,j,k) = CS%land_val(m)
          else
            local_tr(i,j,k) = CS%tr(i,j,k,m)
          endif
        enddo ; enddo ; enddo
      else
        do k=1,nz ; do j=js,je ; do i=is,ie
          local_tr(i,j,k) = CS%tr(i,j,k,m)
        enddo ; enddo ; enddo
      endif ! CS%mask_tracers
      call post_data(CS%id_tracer(m),local_tr,CS%diag)
    endif ! CS%id_tracer(m)>0
    if (CS%id_tr_adx(m)>0) &
      call post_data(CS%id_tr_adx(m),CS%tr_adx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_ady(m)>0) &
      call post_data(CS%id_tr_ady(m),CS%tr_ady(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfx(m)>0) &
      call post_data(CS%id_tr_dfx(m),CS%tr_dfx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfy(m)>0) &
      call post_data(CS%id_tr_dfy(m),CS%tr_dfy(m)%p(:,:,:),CS%diag)
  enddo
  deallocate(local_tr)

end subroutine boundary_impulse_tracer_column_physics

!> Calculate total inventory of tracer
function boundary_impulse_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G
  type(verticalGrid_type),            intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
  real, dimension(:),                 intent(out)   :: stocks
  type(boundary_impulse_tracer_CS),                pointer       :: CS
  character(len=*), dimension(:),     intent(out)   :: names
  character(len=*), dimension(:),     intent(out)   :: units
  integer, optional,                  intent(in)    :: stock_index
  integer                                           :: boundary_impulse_stock
! This function calculates the mass-weighted integral of all tracer stocks,
! returning the number of stocks it has calculated.  If the stock_index
! is present, only the stock corresponding to that coded index is returned.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (out)     stocks - the mass-weighted integrated amount of each tracer,
!                     in kg times concentration units.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_boundary_impulse_tracer.
!  (out)     names - the names of the stocks calculated.
!  (out)     units - the units of the stocks calculated.
!  (in,opt)  stock_index - the coded index of a specific stock being sought.
! Return value: the number of stocks calculated here.

  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  boundary_impulse_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,1
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="boundary_impulse_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      stocks(m) = stocks(m) + CS%tr(i,j,k,m) * &
                           (G%mask2dT(i,j) * G%areaT(i,j) * h(i,j,k))
    enddo ; enddo ; enddo
    stocks(m) = GV%H_to_kg_m2 * stocks(m)
  enddo

  boundary_impulse_stock = CS%ntr

end function boundary_impulse_stock

!> Called if returned if coupler needs to know about tracer, currently unused
subroutine boundary_impulse_tracer_surface_state(state, h, G, CS)
  type(ocean_grid_type),                 intent(in)    :: G
  type(surface),                         intent(inout) :: state
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h
  type(boundary_impulse_tracer_CS),         pointer       :: CS
!   This particular tracer package does not report anything back to the coupler.
! The code that is here is just a rough guide for packages that would.
! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_boundary_impulse_tracer.
  integer :: m, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,CS%ntr
      !   This call loads the surface vlues into the appropriate array in the
      ! coupler-type structure.
      call set_coupler_values(CS%tr(:,:,1,m), state%tr_fields, CS%ind_tr(m), &
                              ind_csurf, is, ie, js, je)
    enddo
  endif

end subroutine boundary_impulse_tracer_surface_state

! Performs finalization of boundary impulse tracer
subroutine boundary_impulse_tracer_end(CS)
  type(boundary_impulse_tracer_CS), pointer :: CS
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    do m=1,CS%ntr
      if (associated(CS%tr_adx(m)%p)) deallocate(CS%tr_adx(m)%p)
      if (associated(CS%tr_ady(m)%p)) deallocate(CS%tr_ady(m)%p)
      if (associated(CS%tr_dfx(m)%p)) deallocate(CS%tr_dfx(m)%p)
      if (associated(CS%tr_dfy(m)%p)) deallocate(CS%tr_dfy(m)%p)
    enddo

    deallocate(CS)
  endif
end subroutine boundary_impulse_tracer_end

!> \namespace boundary_impulse_tracer
!! \section Boundary Impulse Response Tracer and Transit Time Distributions
!! Transit time distributions (TTD) are the Green's function solution of the passive tracer equation between
!! the oceanic surface and interior. The name derives from the idea that the 'age' (e.g. time since last
!! contact with the atmosphere) of a water parcel is best characterized as a distribution of ages
!! because water parcels leaving the surface arrive at a particular interior point at different times.
!! The more commonly used ideal age tracer is the first moment of the TTD, equivalently referred to as the
!! mean age.
!!
!! A boundary impulse response (BIR) is a passive tracer whose surface boundary condition is a rectangle
!! function with width $\Delta t$. In the case of unsteady flow, multiple BIRs, initiated at different
!! times in the model can be used to infer the transit time distribution or Green's function between
!! the oceanic surface and interior. In the case of steady or cyclostationary flow, a single BIR is
!! sufficient.
!!
!! In the References section, both the theoretical discussion of TTDs and BIRs are listed along with
!! modeling studies which have this used framework in scientific investigations
!!
!! \section Run-time parameters
!! -DO_BOUNDARY_IMPULSE_TRACER: Enables the boundary impulse tracer model
!! -IMPULSE_SOURCE_TIME: Length of time that the surface layer acts as a source of the BIR tracer
!!
!! \section References
!! \subsection TTD and BIR Theory
!!  -Holzer, M., and T.M. Hall, 2000: Transit-time and tracer-age distributions in geophysical flows.
!!    J. Atmos. Sci., 57, 3539-3558, doi:10.1175/1520-0469(2000)057<3539:TTATAD>2.0.CO;2.
!!  -T.W.N. Haine, H. Zhang, D.W. Waugh, M. Holzer, On transit-time distributions in unsteady circulation
!!    models, Ocean Modelling, Volume 21, Issues 1–2, 2008, Pages 35-45, ISSN 1463-5003
!!    http://dx.doi.org/10.1016/j.ocemod.2007.11.004.
!! \subsection BIR Modelling applications
!!  -Peacock, S., and M. Maltrud (2006), Transit-time distributions in a global ocean model,
!!    J. Phys. Oceanogr., 36(3), 474–495, doi:10.1175/JPO2860.1.
!!  -Maltrud, M., Bryan, F. & Peacock, Boundary impulse response functions in a century-long eddying global
!!  ocean simulation, S. Environ Fluid Mech (2010) 10: 275. doi:10.1007/s10652-009-9154-3
!!
end module boundary_impulse_tracer
