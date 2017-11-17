!> This module contains the tracer_registry_type and the subroutines
!! that handle registration of tracers and related subroutines.
!! The primary subroutine, register_tracer, is called to indicate the
!! tracers advected and diffused.
module MOM_tracer_registry

! This file is part of MOM6. See LICENSE.md for the license.

! use MOM_diag_mediator, only : diag_ctrl
use MOM_coms,          only : reproducing_sum
use MOM_debugging,     only : hchksum
use MOM_diag_mediator, only : diag_ctrl, register_diag_field, post_data, safe_alloc_ptr
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
public MOM_tracer_chksum, MOM_tracer_chkinv
public register_tracer_diagnostics, post_tracer_diagnostics
public add_tracer_diagnostics, add_tracer_OBC_values
public tracer_registry_init, lock_tracer_registry, tracer_registry_end

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

  character(len=32)               :: name                     !< tracer name used for diagnostics and error messages
  type(vardesc), pointer          :: vd             => NULL() !< metadata describing the tracer
  logical                         :: registry_diags = .false. !< If true, use the registry to set up the
                                                              !! diagnostics associated with this tracer.
  character(len=32)               :: flux_nameroot = ""       !< Short tracer name snippet used construct the
                                                              !! names of flux diagnostics.
  character(len=64)               :: flux_longname = ""       !< A word or phrase used construct the long
                                                              !! names of flux diagnostics.
  real                            :: flux_conversion = 1.0    !< A scaling factor used to convert the fluxes
                                                              !! of this tracer to its desired units.
  character(len=48)               :: flux_units = ""          !< The units for fluxes of this variable.
  integer :: id_tr = -1
  integer :: id_adx = -1, id_ady = -1, id_dfx = -1, id_dfy = -1
  integer :: id_adx_2d = -1, id_ady_2d = -1, id_dfx_2d = -1, id_dfy_2d = -1
  integer :: id_adv_xy = -1, id_adv_xy_2d = -1
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
                           ad_2d_x, ad_2d_y, df_2d_x, df_2d_y, advection_xy, registry_diags, &
                           flux_nameroot, flux_longname, flux_units, flux_conversion)
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
  logical,              optional, intent(in)    :: registry_diags !< If present and true, use the registry for
                                                                !! the diagnostics of this tracer.
  character(len=*),     optional, intent(in)    :: flux_nameroot !< Short tracer name snippet used construct the
                                                                !! names of flux diagnostics.
  character(len=*),     optional, intent(in)    :: flux_longname !< A word or phrase used construct the long
                                                                !! names of flux diagnostics.
  character(len=*),     optional, intent(in)    :: flux_units   !< A scaling factor used to convert the fluxes
                                                                !! of this tracer to its desired units.
  real,                 optional, intent(in)    :: flux_conversion !< A scaling factor used to convert the fluxes
                                                                !! of this tracer to its desired units.
  integer :: ntr
  type(tracer_type) :: temp
  character(len=72) :: longname ! The long name of a variable.
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

  call query_vardesc(Reg%Tr(ntr)%vd, name=Reg%Tr(ntr)%name, longname=longname)

  if (Reg%locked) call MOM_error(FATAL, &
      "MOM register_tracer was called for variable "//trim(Reg%Tr(ntr)%name)//&
      " with a locked tracer registry.")

  Reg%Tr(ntr)%flux_nameroot = Reg%Tr(ntr)%name
  if (present(flux_nameroot)) then
    if (len_trim(flux_nameroot) > 0) Reg%Tr(ntr)%flux_nameroot = flux_nameroot
  endif

  Reg%Tr(ntr)%flux_longname = longname
  if (present(flux_longname)) then
    if (len_trim(flux_longname) > 0) Reg%Tr(ntr)%flux_longname = flux_longname
  endif

  Reg%Tr(ntr)%flux_units = ""
  if (present(flux_units)) Reg%Tr(ntr)%flux_units = flux_units

  Reg%Tr(ntr)%flux_conversion = 1.0
  if (present(flux_conversion)) Reg%Tr(ntr)%flux_conversion = flux_conversion

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

  if (present(registry_diags)) then
    Reg%Tr(ntr)%registry_diags = registry_diags
  endif

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

!> register_tracer_diagnostics does a set of register_diag_field calls for any previously
!! registered in a tracer registry with a value of registry_diags set to .true.
subroutine register_tracer_diagnostics(Reg, Time, diag, G, GV)
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  type(time_type),            intent(in) :: Time !< current model time
  type(diag_ctrl),            intent(in) :: diag !< structure to regulate diagnostic output
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure

  character(len=24) :: name     ! A variable's name in a NetCDF file.
  character(len=24) :: shortnm  ! A shortened version of a variable's name for
                                ! creating additional diagnostics.
  character(len=72) :: longname ! The long name of that tracer variable.
  character(len=72) :: flux_longname ! The tracer name in the long names of fluxes.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for fluxes, either
                                ! [units] m3 s-1 or [units] kg s-1.
  character(len=72) :: cmorname ! The CMOR name of that variable.
  type(tracer_type), pointer :: Tr=>NULL()
  integer :: m, form
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  form = 1

  if (.not. associated(Reg)) call MOM_error(FATAL, "add_tracer_diagnostics: "// &
       "register_tracer must be called before add_tracer_diagnostics")

  do m=1,Reg%ntr ; if (Reg%Tr(m)%registry_diags) then
    Tr => Reg%Tr(m)
    call query_vardesc(Tr%vd, name, units=units, longname=longname, &
                       cmor_field_name=cmorname, caller="register_tracer_diagnostics")
    shortnm = Tr%flux_nameroot
    flux_longname = Tr%flux_longname

    if (len_trim(Tr%flux_units) > 0) then ; flux_units = Tr%flux_units
    elseif (GV%Boussinesq) then ; flux_units = trim(units)//" m3 s-1"
    else ; flux_units = trim(units)//" kg s-1" ; endif

    if (len_trim(cmorname) == 0) then
      Tr%id_tr = register_diag_field("ocean_model", trim(name), diag%axesTL, &
        Time, trim(longname), trim(units))
    else
      Tr%id_tr = register_diag_field("ocean_model", trim(name), diag%axesTL, &
        Time, trim(longname), trim(units), cmor_field_name=cmorname)
    endif
    if (form == 1) then
      Tr%id_adx = register_diag_field("ocean_model", trim(shortnm)//"_adx", &
          diag%axesCuL, Time, trim(flux_longname)//" advective zonal flux" , &
          trim(flux_units))
      Tr%id_ady = register_diag_field("ocean_model", trim(shortnm)//"_ady", &
          diag%axesCvL, Time, trim(flux_longname)//" advective meridional flux" , &
          trim(flux_units))
      Tr%id_dfx = register_diag_field("ocean_model", trim(shortnm)//"_dfx", &
          diag%axesCuL, Time, trim(flux_longname)//" diffusive zonal flux" , &
          trim(flux_units))
      Tr%id_dfy = register_diag_field("ocean_model", trim(shortnm)//"_dfy", &
          diag%axesCvL, Time, trim(flux_longname)//" diffusive zonal flux" , &
          trim(flux_units))
    else
      Tr%id_adx = register_diag_field("ocean_model", trim(shortnm)//"_adx", &
          diag%axesCuL, Time, "Advective (by residual mean) Zonal Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=Tr%flux_conversion)
      Tr%id_ady = register_diag_field("ocean_model", trim(shortnm)//"_ady", &
          diag%axesCvL, Time, "Advective (by residual mean) Meridional Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=Tr%flux_conversion)
      Tr%id_dfx = register_diag_field("ocean_model", trim(shortnm)//"_diffx", &
          diag%axesCuL, Time, "Diffusive Zonal Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=Tr%flux_conversion)
      Tr%id_dfy = register_diag_field("ocean_model", trim(shortnm)//"_diffy", &
          diag%axesCvL, Time, "Diffusive Meridional Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=Tr%flux_conversion)
    endif
    if (Tr%id_adx > 0) call safe_alloc_ptr(Tr%ad_x,IsdB,IedB,jsd,jed,nz)
    if (Tr%id_ady > 0) call safe_alloc_ptr(Tr%ad_y,isd,ied,JsdB,JedB,nz)
    if (Tr%id_dfx > 0) call safe_alloc_ptr(Tr%df_x,IsdB,IedB,jsd,jed,nz)
    if (Tr%id_dfy > 0) call safe_alloc_ptr(Tr%df_y,isd,ied,JsdB,JedB,nz)

    Tr%id_adx_2d = register_diag_field("ocean_model", trim(shortnm)//"_adx_2d", &
        diag%axesCu1, Time, &
        "Vertically Integrated Advective Zonal Flux of "//trim(flux_longname), &
        flux_units, v_extensive=.true., conversion=Tr%flux_conversion)
    Tr%id_ady_2d = register_diag_field("ocean_model", trim(shortnm)//"_ady_2d", &
        diag%axesCv1, Time, &
        "Vertically Integrated Advective Meridional Flux of "//trim(flux_longname), &
        flux_units, v_extensive=.true., conversion=Tr%flux_conversion)
    Tr%id_dfx_2d = register_diag_field("ocean_model", trim(shortnm)//"_diffx_2d", &
        diag%axesCu1, Time, &
        "Vertically Integrated Diffusive Zonal Flux of "//trim(flux_longname), &
        flux_units, v_extensive=.true., conversion=Tr%flux_conversion)
    Tr%id_dfy_2d = register_diag_field("ocean_model", trim(shortnm)//"_diffy_2d", &
        diag%axesCv1, Time, &
        "Vertically Integrated Diffusive Meridional Flux of "//trim(flux_longname), &
        flux_units, v_extensive=.true., conversion=Tr%flux_conversion)

    if (Tr%id_adx_2d > 0) call safe_alloc_ptr(Tr%ad2d_x,IsdB,IedB,jsd,jed)
    if (Tr%id_ady_2d > 0) call safe_alloc_ptr(Tr%ad2d_y,isd,ied,JsdB,JedB)
    if (Tr%id_dfx_2d > 0) call safe_alloc_ptr(Tr%df2d_x,IsdB,IedB,jsd,jed)
    if (Tr%id_dfy_2d > 0) call safe_alloc_ptr(Tr%df2d_y,isd,ied,JsdB,JedB)

    Tr%id_adv_xy = register_diag_field('ocean_model', trim(shortnm)//"_advection_xy", &
        diag%axesTL, Time, &
        'Horizontal convergence of residual mean advective fluxes of '//trim(flux_longname), &
        flux_units, v_extensive=.true., conversion=Tr%flux_conversion)
    Tr%id_adv_xy_2d = register_diag_field('ocean_model', trim(shortnm)//"_advection_xy_2d", &
        diag%axesT1, Time, &
        'Vertical sum of horizontal convergence of residual mean advective fluxes of'//trim(flux_longname), &
        flux_units, conversion=Tr%flux_conversion)
    if ((Tr%id_adv_xy > 0) .or. (Tr%id_adv_xy_2d > 0)) &
      call safe_alloc_ptr(Tr%advection_xy,isd,ied,jsd,jed,nz)

!    call register_Z_tracer(Tr%t, name, longname, units, &
!                           Time, G, diag_to_Z_CSp)
  endif ; enddo

end subroutine register_tracer_diagnostics

!> post_tracer_diagnostics does post_data calls for any diagnostics that are
!! being handled via the tracer registry.
subroutine post_tracer_diagnostics(Reg, diag, G, GV)
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  type(diag_ctrl),            intent(in) :: diag !< structure to regulate diagnostic output
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure

  real    :: work2d(SZI_(G),SZJ_(G))
  type(tracer_type), pointer :: Tr=>NULL()
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  do m=1,Reg%ntr ; if (Reg%Tr(m)%registry_diags) then
    Tr => Reg%Tr(m)
    if (Tr%id_tr > 0) call post_data(Tr%id_tr, Tr%t, diag)
    if (Tr%id_adx > 0) call post_data(Tr%id_adx, Tr%ad_x, diag)
    if (Tr%id_ady > 0) call post_data(Tr%id_ady, Tr%ad_y, diag)
    if (Tr%id_dfx > 0) call post_data(Tr%id_dfx, Tr%df_x, diag)
    if (Tr%id_dfy > 0) call post_data(Tr%id_dfy, Tr%df_y, diag)
    if (Tr%id_adx_2d > 0) call post_data(Tr%id_adx_2d, Tr%ad2d_x, diag)
    if (Tr%id_ady_2d > 0) call post_data(Tr%id_ady_2d, Tr%ad2d_y, diag)
    if (Tr%id_dfx_2d > 0) call post_data(Tr%id_dfx_2d, Tr%df2d_x, diag)
    if (Tr%id_dfy_2d > 0) call post_data(Tr%id_dfy_2d, Tr%df2d_y, diag)
    if (Tr%id_adv_xy > 0) call post_data(Tr%id_adv_xy, Tr%advection_xy, diag)
    if (Tr%id_adv_xy_2d > 0) then
      work2d(:,:) = 0.0
      do k=1,nz ; do j=js,je ; do i=is,ie
        work2d(i,j) = work2d(i,j) + Tr%advection_xy(i,j,k)
      enddo ; enddo ; enddo
      call post_data(Tr%id_adv_xy_2d, Tr%advection_xy, diag)
    endif
  endif ; enddo
end subroutine post_tracer_diagnostics

!> This subroutine writes out chksums for tracers.
subroutine MOM_tracer_chksum(mesg, Tr, ntr, G)
  character(len=*),         intent(in) :: mesg   !< message that appears on the chksum lines
  type(tracer_type),        intent(in) :: Tr(:)  !< array of all of registered tracers
  integer,                  intent(in) :: ntr    !< number of registered tracers
  type(ocean_grid_type),    intent(in) :: G      !< ocean grid structure

  integer :: is, ie, js, je, nz
  integer :: i, j, k, m

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  do m=1,ntr
    call hchksum(Tr(m)%t, mesg//trim(Tr(m)%name), G%HI)
  enddo

end subroutine MOM_tracer_chksum

!> Calculates and prints the global inventory of all tracers in the registry.
subroutine MOM_tracer_chkinv(mesg, G, h, Tr, ntr)
  character(len=*),                         intent(in) :: mesg   !< message that appears on the chksum lines
  type(ocean_grid_type),                    intent(in) :: G      !< ocean grid structure
  type(tracer_type),                        intent(in) :: Tr(:)  !< array of all of registered tracers
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h      !< Layer thicknesses
  integer,                                  intent(in) :: ntr    !< number of registered tracers

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: tr_inv !< Tracer inventory
  real :: total_inv
  integer :: is, ie, js, je, nz
  integer :: i, j, k, m

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  do m=1,ntr
    do k=1,nz ; do j=js,je ; do i=is,ie
      tr_inv(i,j,k) = Tr(m)%t(i,j,k)*h(i,j,k)*G%areaT(i,j)*G%mask2dT(i,j)
    enddo ; enddo ; enddo
    total_inv = reproducing_sum(tr_inv, is, ie, js, je)
    if (is_root_pe()) write(0,'(A,1X,A5,1X,ES25.16,1X,A)') "h-point: inventory", Tr(m)%name, total_inv, mesg
  enddo

end subroutine MOM_tracer_chkinv

!> Initialize the tracer registry.
subroutine tracer_registry_init(param_file, Reg)
  type(param_file_type),      intent(in) :: param_file !< open file to parse for model parameters
  type(tracer_registry_type), pointer    :: Reg        !< pointer to tracer registry

  integer, save :: init_calls = 0

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_tracer_registry" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (.not.associated(Reg)) then ; allocate(Reg)
  else ; return ; endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

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
