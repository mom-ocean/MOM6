!> This module contains the tracer_registry_type and the subroutines
!! that handle registration of tracers and related subroutines.
!! The primary subroutine, register_tracer, is called to indicate the
!! tracers advected and diffused.
module MOM_tracer_registry

! This file is part of MOM6. See LICENSE.md for the license.

! use MOM_diag_mediator, only : diag_ctrl
use MOM_debugging,     only : hchksum
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_hor_index,     only : hor_index_type
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : vardesc, query_vardesc
use MOM_time_manager,  only : time_type
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_tracer
public tracer_registry_init
public MOM_tracer_chksum
public add_tracer_diagnostics
public add_tracer_OBC_values
public lock_tracer_registry
public tracer_registry_end

!> The tracer type
type, public :: tracer_type

  real, dimension(:,:,:), pointer :: t              => NULL() !< tracer concentration array
  real                            :: OBC_inflow_conc=  0.0    !< tracer concentration for generic inflows
  real, dimension(:,:,:), pointer :: OBC_in_u       => NULL() !< structured values for flow into the domain
                                                              !! specified in OBCs through u-face of cell
  real, dimension(:,:,:), pointer :: OBC_in_v       => NULL() !< structured values for flow into the domain
                                                              !! specified in OBCs through v-face of cell

  real, dimension(:,:,:), pointer :: ad_x           => NULL() !< diagnostic array for x-advective tracer flux
  real, dimension(:,:,:), pointer :: ad_y           => NULL() !< diagnostic array for y-advective tracer flux
  real, dimension(:,:),   pointer :: ad2d_x         => NULL() !< diagnostic vertical sum x-advective tracer flux
                                                              !! in units of (conc * m3/s or conc * kg/s)
  real, dimension(:,:),   pointer :: ad2d_y         => NULL() !< diagnostic vertical sum y-advective tracer flux
                                                              !! in units of (conc * m3/s or conc * kg/s)

  real, dimension(:,:,:), pointer :: df_x           => NULL() !< diagnostic array for x-diffusive tracer flux
  real, dimension(:,:,:), pointer :: df_y           => NULL() !< diagnostic array for y-diffusive tracer flux
  real, dimension(:,:),   pointer :: df2d_x         => NULL() !< diagnostic vertical sum x-diffusive flux
                                                              !! in units of (conc * m3/s or conc * kg/s)
  real, dimension(:,:),   pointer :: df2d_y         => NULL() !< diagnostic vertical sum y-diffusive flux
                                                              !! in units of (conc * m3/s or conc * kg/s)

  real, dimension(:,:,:), pointer :: advection_xy   => NULL() !< convergence of lateral advective tracer fluxes

  character(len=32)               :: name                     !< tracer name used for error messages
  type(vardesc), pointer          :: vd             => NULL() !< metadata describing the tracer

end type tracer_type

!> Type to carry basic tracer information
type, public :: tracer_registry_type
  integer                  :: ntr = 0           !< number of registered tracers
  type(tracer_type)        :: Tr(MAX_FIELDS_)   !< array of registered tracers
! type(diag_ctrl), pointer :: diag              !< structure to regulate timing of diagnostics
  logical                  :: locked = .false.  !< New tracers may be registered if locked=.false.
                                                !! When locked=.true.,no more tracers can be registered,
                                                !! at which point common diagnostics can be set up
                                                !! for the registered tracers.
end type tracer_registry_type

contains

!> This subroutine registers a tracer to be advected and laterally diffused.
subroutine register_tracer(tr1, tr_desc, param_file, HI, GV, Reg, tr_desc_ptr, ad_x, ad_y,&
                           df_x, df_y, OBC_inflow, OBC_in_u, OBC_in_v,            &
                           ad_2d_x, ad_2d_y, df_2d_x, df_2d_y, advection_xy)
  type(hor_index_type),           intent(in)    :: HI           !< horizontal index type
  type(verticalGrid_type),        intent(in)    :: GV           !< ocean vertical grid structure
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), target :: tr1    !< pointer to the tracer (concentration units)
  type(vardesc),         intent(in)             :: tr_desc      !< metadata about the tracer
  type(param_file_type), intent(in)             :: param_file   !< file to parse for  model parameter values
  type(tracer_registry_type), pointer           :: Reg          !< pointer to the tracer registry
  type(vardesc), target, optional               :: tr_desc_ptr  !< A target that can be used to set a pointer to the
                                                                !! stored value of tr%tr_desc.  This target must be an
                                                                !! enduring part of the control structure, because the tracer
                                                                !! registry will use this memory, but it also means that any
                                                                !! updates to this structure in the calling module will be
                                                                !! available subsequently to the tracer registry.
  real, pointer, dimension(:,:,:), optional     :: ad_x         !< diagnostic x-advective flux (CONC m3/s or CONC*kg/s)
  real, pointer, dimension(:,:,:), optional     :: ad_y         !< diagnostic y-advective flux (CONC m3/s or CONC*kg/s)
  real, pointer, dimension(:,:,:), optional     :: df_x         !< diagnostic x-diffusive flux (CONC m3/s or CONC*kg/s)
  real, pointer, dimension(:,:,:), optional     :: df_y         !< diagnostic y-diffusive flux (CONC m3/s or CONC*kg/s)

  real, intent(in),                optional     :: OBC_inflow   !< the tracer for all inflows via OBC for which OBC_in_u
                                                                !! or OBC_in_v are not specified (units of tracer CONC)
  real, pointer, dimension(:,:,:), optional     :: OBC_in_u     !< tracer at inflows through u-faces of
                                                                !! tracer cells (units of tracer CONC)
  real, pointer, dimension(:,:,:), optional     :: OBC_in_v     !< tracer at inflows through v-faces of
                                                                !! tracer cells (units of tracer CONC)

  real, dimension(:,:),   pointer, optional     :: ad_2d_x      !< vert sum of diagnostic x-advect flux (CONC m3/s or CONC*kg/s)
  real, dimension(:,:),   pointer, optional     :: ad_2d_y      !< vert sum of diagnostic y-advect flux (CONC m3/s or CONC*kg/s)
  real, dimension(:,:),   pointer, optional     :: df_2d_x      !< vert sum of diagnostic x-diffuse flux (CONC m3/s or CONC*kg/s)
  real, dimension(:,:),   pointer, optional     :: df_2d_y      !< vert sum of diagnostic y-diffuse flux (CONC m3/s or CONC*kg/s)

  real, pointer, dimension(:,:,:), optional     :: advection_xy !< convergence of lateral advective tracer fluxes

  integer :: ntr
  type(tracer_type) :: temp
  character(len=256) :: mesg    ! Message for error messages.

  if (.not. associated(Reg)) call tracer_registry_init(param_file, Reg)

  if (Reg%ntr>=MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ in MOM_memory.h to at least ",I3," to allow for &
        &all the tracers being registered via register_tracer.")') Reg%ntr+1
    call MOM_error(FATAL,"MOM register_tracer: "//mesg)
  endif
  Reg%ntr = Reg%ntr + 1
  ntr     = Reg%ntr

  if (present(tr_desc_ptr)) then
    Reg%Tr(ntr)%vd => tr_desc_ptr
  else
    allocate(Reg%Tr(ntr)%vd) ; Reg%Tr(ntr)%vd = tr_desc
  endif

  call query_vardesc(Reg%Tr(ntr)%vd, name=Reg%Tr(ntr)%name)

  if (Reg%locked) call MOM_error(FATAL, &
      "MOM register_tracer was called for variable "//trim(Reg%Tr(ntr)%name)//&
      " with a locked tracer registry.")

  Reg%Tr(ntr)%t => tr1

  if (present(ad_x)) then ; if (associated(ad_x)) Reg%Tr(ntr)%ad_x => ad_x ; endif
  if (present(ad_y)) then ; if (associated(ad_y)) Reg%Tr(ntr)%ad_y => ad_y ; endif
  if (present(df_x)) then ; if (associated(df_x)) Reg%Tr(ntr)%df_x => df_x ; endif
  if (present(df_y)) then ; if (associated(df_y)) Reg%Tr(ntr)%df_y => df_y ; endif
  if (present(OBC_inflow)) Reg%Tr(ntr)%OBC_inflow_conc = OBC_inflow
  if (present(OBC_in_u)) then ; if (associated(OBC_in_u)) &
                                    Reg%Tr(ntr)%OBC_in_u => OBC_in_u ; endif
  if (present(OBC_in_v)) then ; if (associated(OBC_in_v)) &
                                    Reg%Tr(ntr)%OBC_in_v => OBC_in_v ; endif
  if (present(ad_2d_x)) then ; if (associated(ad_2d_x)) Reg%Tr(ntr)%ad2d_x => ad_2d_x ; endif
  if (present(ad_2d_y)) then ; if (associated(ad_2d_y)) Reg%Tr(ntr)%ad2d_y => ad_2d_y ; endif
  if (present(df_2d_x)) then ; if (associated(df_2d_x)) Reg%Tr(ntr)%df2d_x => df_2d_x ; endif

  if (present(advection_xy)) then ; if (associated(advection_xy)) Reg%Tr(ntr)%advection_xy => advection_xy ; endif

end subroutine register_tracer


!> This subroutine locks the tracer registry to prevent the addition of more
!! tracers.  After locked=.true., can then register common diagnostics.
subroutine lock_tracer_registry(Reg)
  type(tracer_registry_type), pointer    :: Reg    !< pointer to the tracer registry

  if (.not. associated(Reg)) call MOM_error(WARNING, &
    "lock_tracer_registry called with an unassocaited registry.")

  Reg%locked = .True.

end subroutine lock_tracer_registry


!> This subroutine adds open boundary condition concentrations for a tracer that
!! has previously been registered by a call to register_tracer.
subroutine add_tracer_OBC_values(name, Reg, OBC_inflow, OBC_in_u, OBC_in_v)
  character(len=*), intent(in)               :: name        !< tracer name for which the diagnostic points
  type(tracer_registry_type), pointer        :: Reg         !< pointer to the tracer registry
  real, intent(in), optional                 :: OBC_inflow  !< tracer value for all inflows via the OBC
                                                            !! for which OBC_in_u or OBC_in_v are
                                                            !! not specified (same units as tracer CONC)
  real, pointer, dimension(:,:,:), optional  :: OBC_in_u    !< tracer at inflows through u-face of tracer cells
                                                            !! (same units as tracer CONC)
  real, pointer, dimension(:,:,:), optional  :: OBC_in_v    !< tracer at inflows through v-face of tracer cells
                                                            !! (same units as tracer CONC)

  integer :: m

  if (.not. associated(Reg)) call MOM_error(FATAL, "add_tracer_OBC_values :"// &
       "register_tracer must be called before add_tracer_OBC_values")

  do m=1,Reg%ntr ; if (Reg%Tr(m)%name == trim(name)) exit ; enddo

  if (m <= Reg%ntr) then
    if (present(OBC_inflow)) Reg%Tr(m)%OBC_inflow_conc = OBC_inflow
    if (present(OBC_in_u)) then ; if (associated(OBC_in_u)) &
                                      Reg%Tr(m)%OBC_in_u => OBC_in_u ; endif
    if (present(OBC_in_v)) then ; if (associated(OBC_in_v)) &
                                      Reg%Tr(m)%OBC_in_v => OBC_in_v ; endif
  else
    call MOM_error(FATAL, "MOM_tracer: register_tracer must be called for "//&
             trim(name)//" before add_tracer_OBC_values is called for it.")
  endif

end subroutine add_tracer_OBC_values


!> This subroutine adds diagnostic arrays for a tracer that has
!! previously been registered by a call to register_tracer.
subroutine add_tracer_diagnostics(name, Reg, ad_x, ad_y, df_x, df_y, &
                                  ad_2d_x, ad_2d_y, df_2d_x, df_2d_y,&
                                  advection_xy)
  character(len=*), intent(in)              :: name         !< name of the tracer for which the diagnostic points
  type(tracer_registry_type), pointer       :: Reg          !< pointer to the tracer registry
  real, dimension(:,:,:), pointer, optional :: ad_x         !< diagnostic x-advective flux (CONC m3/s or CONC*kg/s)
  real, dimension(:,:,:), pointer, optional :: ad_y         !< diagnostic y-advective flux (CONC m3/s or CONC*kg/s)
  real, dimension(:,:,:), pointer, optional :: df_x         !< diagnostic x-diffusive flux (CONC m3/s or CONC*kg/s)
  real, dimension(:,:,:), pointer, optional :: df_y         !< diagnostic y-diffusive flux (CONC m3/s or CONC*kg/s)
  real, dimension(:,:),   pointer, optional :: ad_2d_x      !< vert sum of diagnostic x-advect flux (CONC m3/s or CONC*kg/s)
  real, dimension(:,:),   pointer, optional :: ad_2d_y      !< vert sum of diagnostic y-advect flux (CONC m3/s or CONC*kg/s)
  real, dimension(:,:),   pointer, optional :: df_2d_x      !< vert sum of diagnostic x-diffuse flux (CONC m3/s or CONC*kg/s)
  real, dimension(:,:),   pointer, optional :: df_2d_y      !< vert sum of diagnostic y-diffuse flux (CONC m3/s or CONC*kg/s)

  real, dimension(:,:,:), pointer, optional :: advection_xy !< convergence of lateral advective tracer fluxes

  integer :: m

  if (.not. associated(Reg)) call MOM_error(FATAL, "add_tracer_diagnostics: "// &
       "register_tracer must be called before add_tracer_diagnostics")

  do m=1,Reg%ntr ; if (Reg%Tr(m)%name == trim(name)) exit ; enddo

  if (m <= Reg%ntr) then
    if (present(ad_x)) then ; if (associated(ad_x)) Reg%Tr(m)%ad_x => ad_x ; endif
    if (present(ad_y)) then ; if (associated(ad_y)) Reg%Tr(m)%ad_y => ad_y ; endif
    if (present(df_x)) then ; if (associated(df_x)) Reg%Tr(m)%df_x => df_x ; endif
    if (present(df_y)) then ; if (associated(df_y)) Reg%Tr(m)%df_y => df_y ; endif

    if (present(ad_2d_x)) then ; if (associated(ad_2d_x)) Reg%Tr(m)%ad2d_x => ad_2d_x ; endif
    if (present(ad_2d_y)) then ; if (associated(ad_2d_y)) Reg%Tr(m)%ad2d_y => ad_2d_y ; endif
    if (present(df_2d_x)) then ; if (associated(df_2d_x)) Reg%Tr(m)%df2d_x => df_2d_x ; endif
    if (present(df_2d_y)) then ; if (associated(df_2d_y)) Reg%Tr(m)%df2d_y => df_2d_y ; endif

    if (present(advection_xy)) then ; if (associated(advection_xy)) Reg%Tr(m)%advection_xy => advection_xy ; endif

  else

    call MOM_error(FATAL, "MOM_tracer: register_tracer must be called for "//&
             trim(name)//" before add_tracer_diagnostics is called for it.")
  endif

end subroutine add_tracer_diagnostics

!> This subroutine writes out chksums for thermodynamic state variables.
subroutine MOM_tracer_chksum(mesg, Tr, ntr, G)
  character(len=*),         intent(in) :: mesg   !< message that appears on the chksum lines
  type(tracer_type),        intent(in) :: Tr(:)  !< array of all of registered tracers
  integer,                  intent(in) :: ntr    !< number of registered tracers
  type(ocean_grid_type),    intent(in) :: G      !< ocean grid structure

  integer :: is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  do m=1,ntr
    call hchksum(Tr(m)%t, mesg//trim(Tr(m)%name), G%HI)
  enddo

end subroutine MOM_tracer_chksum


!> This routine include declares and sets the variable "version".
subroutine tracer_registry_init(param_file, Reg)
  type(param_file_type),      intent(in) :: param_file !< open file to parse for model parameters
  type(tracer_registry_type), pointer    :: Reg        !< pointer to tracer registry

  integer, save :: init_calls = 0

#include "version_variable.h"
  character(len=40)  :: mod = "MOM_tracer_registry" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (.not.associated(Reg)) then ; allocate(Reg)
  else ; return ; endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")

  init_calls = init_calls + 1
  if (init_calls > 1) then
    write(mesg,'("tracer_registry_init called ",I3, &
      &" times with different registry pointers.")') init_calls
    if (is_root_pe()) call MOM_error(WARNING,"MOM_tracer"//mesg)
  endif

end subroutine tracer_registry_init


!> This routine closes the tracer registry module.
subroutine tracer_registry_end(Reg)
  type(tracer_registry_type), pointer :: Reg
  if (associated(Reg)) deallocate(Reg)
end subroutine tracer_registry_end

end module MOM_tracer_registry
