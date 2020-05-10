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
use MOM_diag_mediator, only : diag_grid_storage
use MOM_diag_mediator, only : diag_copy_storage_to_diag, diag_save_grids, diag_restore_grids
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_hor_index,     only : hor_index_type
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : vardesc, query_vardesc, cmor_long_std
use MOM_restart,       only : register_restart_field, MOM_restart_CS
use MOM_string_functions, only : lowercase
use MOM_time_manager,  only : time_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_tracer
public MOM_tracer_chksum, MOM_tracer_chkinv
public register_tracer_diagnostics, post_tracer_diagnostics, post_tracer_transport_diagnostics
public preALE_tracer_diagnostics, postALE_tracer_diagnostics
public tracer_registry_init, lock_tracer_registry, tracer_registry_end
public tracer_name_lookup

!> The tracer type
type, public :: tracer_type

  real, dimension(:,:,:), pointer :: t              => NULL() !< tracer concentration array [conc]
! real                            :: OBC_inflow_conc=  0.0    !< tracer concentration for generic inflows [conc]
! real, dimension(:,:,:), pointer :: OBC_in_u       => NULL() !< structured values for flow into the domain
!                                                             !! specified in OBCs through u-face of cell
! real, dimension(:,:,:), pointer :: OBC_in_v       => NULL() !< structured values for flow into the domain
!                                                             !! specified in OBCs through v-face of cell

  real, dimension(:,:,:), pointer :: ad_x           => NULL() !< diagnostic array for x-advective tracer flux
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:,:), pointer :: ad_y           => NULL() !< diagnostic array for y-advective tracer flux
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   pointer :: ad2d_x         => NULL() !< diagnostic vertical sum x-advective tracer flux
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   pointer :: ad2d_y         => NULL() !< diagnostic vertical sum y-advective tracer flux
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]

  real, dimension(:,:,:), pointer :: df_x           => NULL() !< diagnostic array for x-diffusive tracer flux
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:,:), pointer :: df_y           => NULL() !< diagnostic array for y-diffusive tracer flux
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:,:), pointer :: lbd_dfx       => NULL()  !< diagnostic array for x-diffusive tracer flux
                                                              !! [conc H m2 s-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:,:), pointer :: lbd_dfy       => NULL()  !< diagnostic array for y-diffusive tracer flux
                                                              !! [conc H m2 s-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:), pointer :: lbd_dfx_2d       => NULL() !< diagnostic array for x-diffusive tracer flux
                                                              !! [conc H m2 s-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:), pointer :: lbd_dfy_2d       => NULL() !< diagnostic array for y-diffusive tracer flux
                                                              !! [conc H m2 s-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:), pointer :: lbd_bulk_df_x       => NULL() !< diagnostic array for x-diffusive tracer flux
                                                              !! [conc H m2 s-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:), pointer :: lbd_bulk_df_y       => NULL() !< diagnostic array for y-diffusive tracer flux
                                                              !! [conc H m2 s-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   pointer :: df2d_x         => NULL() !< diagnostic vertical sum x-diffusive flux
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   pointer :: df2d_y         => NULL() !< diagnostic vertical sum y-diffusive flux
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
!  real, dimension(:,:),   pointer :: df2d_conc_x    => NULL() !< diagnostic vertical sum x-diffusive content flux
!                                                              !! [conc H m2 s-1 ~> conc m3 s-1 or conc kg s-1]
!  real, dimension(:,:),   pointer :: df2d_conc_y    => NULL() !< diagnostic vertical sum y-diffusive content flux
!                                                              !! [conc H m2 s-1 ~> conc m3 s-1 or conc kg s-1]

  real, dimension(:,:,:), pointer :: advection_xy   => NULL() !< convergence of lateral advective tracer fluxes
                                                              !! [conc H T-1 ~> conc m s-1 or conc kg m-2 s-1]
!  real, dimension(:,:,:), pointer :: diff_cont_xy   => NULL() !< convergence of lateral diffusive tracer fluxes
!                                                              !! [conc H s-1 ~> conc m s-1 or conc kg m-2 s-1]
!  real, dimension(:,:,:), pointer :: diff_conc_xy   => NULL() !< convergence of lateral diffusive tracer fluxes
!                                                              !! expressed as a change in concentration [conc s-1]
  real, dimension(:,:,:), pointer :: t_prev         => NULL() !< tracer concentration array at a previous
                                                              !! timestep used for diagnostics [conc]
  real, dimension(:,:,:), pointer :: Trxh_prev      => NULL() !< layer integrated tracer concentration array
                                                              !! at a previous timestep used for diagnostics

  character(len=32)               :: name                     !< tracer name used for diagnostics and error messages
  character(len=64)               :: units                    !< Physical dimensions of the tracer concentration
  character(len=240)              :: longname                 !< Long name of the variable
!  type(vardesc), pointer          :: vd             => NULL() !< metadata describing the tracer
  logical                         :: registry_diags = .false. !< If true, use the registry to set up the
                                                              !! diagnostics associated with this tracer.
  character(len=64)               :: cmor_name                !< CMOR name of this tracer
  character(len=64)               :: cmor_units               !< CMOR physical dimensions of the tracer
  character(len=240)              :: cmor_longname            !< CMOR long name of the tracer
  character(len=32)               :: flux_nameroot = ""       !< Short tracer name snippet used construct the
                                                              !! names of flux diagnostics.
  character(len=64)               :: flux_longname = ""       !< A word or phrase used construct the long
                                                              !! names of flux diagnostics.
  real                            :: flux_scale= 1.0          !< A scaling factor used to convert the fluxes
                                                              !! of this tracer to its desired units.
  character(len=48)               :: flux_units = ""          !< The units for fluxes of this variable.
  character(len=48)               :: conv_units = ""          !< The units for the flux convergence of this tracer.
  real                            :: conv_scale = 1.0         !< A scaling factor used to convert the flux
                                                              !! convergence of this tracer to its desired units.
  character(len=48)               :: cmor_tendprefix = ""     !< The CMOR variable prefix for tendencies of this
                                                              !! tracer, required because CMOR does not follow any
                                                              !! discernable pattern for these names.
  integer :: ind_tr_squared = -1 !< The tracer registry index for the square of this tracer

  !### THESE CAPABILITIES HAVE NOT YET BEEN IMPLEMENTED.
  ! logical :: advect_tr = .true.     !< If true, this tracer should be advected
  ! logical :: hordiff_tr = .true.    !< If true, this tracer should experience epineutral diffusion
  logical :: remap_tr = .true.      !< If true, this tracer should be vertically remapped

  integer :: diag_form = 1  !< An integer indicating which template is to be used to label diagnostics.
  !>@{ Diagnostic IDs
  integer :: id_tr = -1
  integer :: id_adx = -1, id_ady = -1, id_dfx = -1, id_dfy = -1
  integer :: id_lbd_bulk_dfx = -1, id_lbd_bulk_dfy = -1, id_lbd_dfx = -1, id_lbd_dfy = -1
  integer :: id_lbd_dfx_2d = -1  , id_lbd_dfy_2d = -1
  integer :: id_adx_2d = -1, id_ady_2d = -1, id_dfx_2d = -1, id_dfy_2d = -1
  integer :: id_adv_xy = -1, id_adv_xy_2d = -1
  integer :: id_dfxy_cont = -1, id_dfxy_cont_2d = -1, id_dfxy_conc = -1
  integer :: id_lbdxy_cont = -1, id_lbdxy_cont_2d = -1, id_lbdxy_conc = -1
  integer :: id_remap_conc = -1, id_remap_cont = -1, id_remap_cont_2d = -1
  integer :: id_tendency = -1, id_trxh_tendency = -1, id_trxh_tendency_2d = -1
  integer :: id_tr_vardec = -1
  !>@}
end type tracer_type

!> Type to carry basic tracer information
type, public :: tracer_registry_type
  integer                  :: ntr = 0           !< number of registered tracers
  type(tracer_type)        :: Tr(MAX_FIELDS_)   !< array of registered tracers
! type(diag_ctrl), pointer :: diag              !< structure to regulate timing of diagnostics
  logical                  :: locked = .false.  !< New tracers may be registered if locked=.false.
                                                !! When locked=.true., no more tracers can be registered,
                                                !! at which point common diagnostics can be set up
                                                !! for the registered tracers
end type tracer_registry_type

contains

!> This subroutine registers a tracer to be advected and laterally diffused.
subroutine register_tracer(tr_ptr, Reg, param_file, HI, GV, name, longname, units, &
                           cmor_name, cmor_units, cmor_longname, tr_desc, &
                           OBC_inflow, OBC_in_u, OBC_in_v, ad_x, ad_y, df_x, df_y, &
                           ad_2d_x, ad_2d_y, df_2d_x, df_2d_y, advection_xy, registry_diags, &
                           flux_nameroot, flux_longname, flux_units, flux_scale, &
                           convergence_units, convergence_scale, cmor_tendprefix, diag_form, &
                           restart_CS, mandatory)
  type(hor_index_type),           intent(in)    :: HI           !< horizontal index type
  type(verticalGrid_type),        intent(in)    :: GV           !< ocean vertical grid structure
  type(tracer_registry_type),     pointer       :: Reg          !< pointer to the tracer registry
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                                  target        :: tr_ptr       !< target or pointer to the tracer array
  type(param_file_type), intent(in)             :: param_file   !< file to parse for model parameter values
  character(len=*),     optional, intent(in)    :: name         !< Short tracer name
  character(len=*),     optional, intent(in)    :: longname     !< The long tracer name
  character(len=*),     optional, intent(in)    :: units        !< The units of this tracer
  character(len=*),     optional, intent(in)    :: cmor_name    !< CMOR name
  character(len=*),     optional, intent(in)    :: cmor_units   !< CMOR physical dimensions of variable
  character(len=*),     optional, intent(in)    :: cmor_longname !< CMOR long name
  type(vardesc),        optional, intent(in)    :: tr_desc      !< A structure with metadata about the tracer

  real,                 optional, intent(in)    :: OBC_inflow   !< the tracer for all inflows via OBC for which OBC_in_u
                                                                !! or OBC_in_v are not specified (units of tracer CONC)
  real, dimension(:,:,:), optional, pointer     :: OBC_in_u     !< tracer at inflows through u-faces of
                                                                !! tracer cells (units of tracer CONC)
  real, dimension(:,:,:), optional, pointer     :: OBC_in_v     !< tracer at inflows through v-faces of
                                                                !! tracer cells (units of tracer CONC)

  ! The following are probably not necessary if registry_diags is present and true.
  real, dimension(:,:,:), optional, pointer     :: ad_x         !< diagnostic x-advective flux
                                                                !! [conc H L2 T-1 ~> CONC m3 s-1 or CONC kg s-1]
  real, dimension(:,:,:), optional, pointer     :: ad_y         !< diagnostic y-advective flux
                                                                !! [conc H L2 T-1 ~> CONC m3 s-1 or CONC kg s-1]
  real, dimension(:,:,:), optional, pointer     :: df_x         !< diagnostic x-diffusive flux
                                                                !! [conc H L2 T-1 ~> CONC m3 s-1 or CONC kg s-1]
  real, dimension(:,:,:), optional, pointer     :: df_y         !< diagnostic y-diffusive flux
                                                                !! [conc H L2 T-1 ~> CONC m3 s-1 or CONC kg s-1]
  real, dimension(:,:),   optional, pointer     :: ad_2d_x      !< vert sum of diagnostic x-advect flux
                                                                !! [conc H L2 T-1 ~> CONC m3 s-1 or CONC kg s-1]
  real, dimension(:,:),   optional, pointer     :: ad_2d_y      !< vert sum of diagnostic y-advect flux
                                                                !! [conc H L2 T-1 ~> CONC m3 s-1 or CONC kg s-1]
  real, dimension(:,:),   optional, pointer     :: df_2d_x      !< vert sum of diagnostic x-diffuse flux
                                                                !! [conc H L2 T-1 ~> CONC m3 s-1 or CONC kg s-1]
  real, dimension(:,:),   optional, pointer     :: df_2d_y      !< vert sum of diagnostic y-diffuse flux
                                                                !! [conc H L2 T-1 ~> CONC m3 s-1 or CONC kg s-1]

  real, dimension(:,:,:), optional, pointer     :: advection_xy !< convergence of lateral advective tracer fluxes
  logical,              optional, intent(in)    :: registry_diags !< If present and true, use the registry for
                                                                !! the diagnostics of this tracer.
  character(len=*),     optional, intent(in)    :: flux_nameroot !< Short tracer name snippet used construct the
                                                                !! names of flux diagnostics.
  character(len=*),     optional, intent(in)    :: flux_longname !< A word or phrase used construct the long
                                                                !! names of flux diagnostics.
  character(len=*),     optional, intent(in)    :: flux_units   !< The units for the fluxes of this tracer.
  real,                 optional, intent(in)    :: flux_scale   !< A scaling factor used to convert the fluxes
                                                                !! of this tracer to its desired units.
  character(len=*),     optional, intent(in)    :: convergence_units !< The units for the flux convergence of
                                                                !! this tracer.
  real,                 optional, intent(in)    :: convergence_scale !< A scaling factor used to convert the flux
                                                                !! convergence of this tracer to its desired units.
  character(len=*),     optional, intent(in)    :: cmor_tendprefix !< The CMOR name for the layer-integrated
                                                                !! tendencies of this tracer.
  integer,              optional, intent(in)    :: diag_form    !< An integer (1 or 2, 1 by default) indicating the
                                                                !! character string template to use in
                                                                !! labeling diagnostics
  type(MOM_restart_CS), optional, pointer       :: restart_CS   !< A pointer to the restart control structure
                                                                !! this tracer will be registered for
                                                                !! restarts if this argument is present
  logical,              optional, intent(in)    :: mandatory    !< If true, this tracer must be read
                                                                !! from a restart file.

  logical :: mand
  type(tracer_type), pointer :: Tr=>NULL()
  character(len=256) :: mesg    ! Message for error messages.

  if (.not. associated(Reg)) call tracer_registry_init(param_file, Reg)

  if (Reg%ntr>=MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ in MOM_memory.h to at least ",I3," to allow for &
        &all the tracers being registered via register_tracer.")') Reg%ntr+1
    call MOM_error(FATAL,"MOM register_tracer: "//mesg)
  endif
  Reg%ntr = Reg%ntr + 1

  Tr => Reg%Tr(Reg%ntr)

  if (present(name)) then
    Tr%name = name
    Tr%longname = name ; if (present(longname)) Tr%longname = longname
    Tr%units = "Conc" ; if (present(units)) Tr%units = units

    Tr%cmor_name = ""
    if (present(cmor_name)) Tr%cmor_name = cmor_name

    Tr%cmor_units = Tr%units
    if (present(cmor_units)) Tr%cmor_units = cmor_units

    Tr%cmor_longname = ""
    if (present(cmor_longname)) Tr%cmor_longname = cmor_longname

    if (present(tr_desc)) call MOM_error(WARNING, "MOM register_tracer: "//&
      "It is a bad idea to use both name and tr_desc when registring "//trim(name))
  elseif (present(tr_desc)) then
    call query_vardesc(tr_desc, name=Tr%name, units=Tr%units, &
                       longname=Tr%longname, cmor_field_name=Tr%cmor_name, &
                       cmor_longname=Tr%cmor_longname, caller="register_tracer")
    Tr%cmor_units = Tr%units
  else
    call MOM_error(FATAL,"MOM register_tracer: Either name or "//&
                   "tr_desc must be present when registering a tracer.")
  endif

  if (Reg%locked) call MOM_error(FATAL, &
      "MOM register_tracer was called for variable "//trim(Tr%name)//&
      " with a locked tracer registry.")

  Tr%flux_nameroot = Tr%name
  if (present(flux_nameroot)) then
    if (len_trim(flux_nameroot) > 0) Tr%flux_nameroot = flux_nameroot
  endif

  Tr%flux_longname = Tr%longname
  if (present(flux_longname)) then
    if (len_trim(flux_longname) > 0) Tr%flux_longname = flux_longname
  endif

  Tr%flux_units = ""
  if (present(flux_units)) Tr%flux_units = flux_units

  Tr%flux_scale = 1.0
  if (present(flux_scale)) Tr%flux_scale = flux_scale

  Tr%conv_units = ""
  if (present(convergence_units)) Tr%conv_units = convergence_units

  Tr%cmor_tendprefix = ""
  if (present(cmor_tendprefix)) Tr%cmor_tendprefix = cmor_tendprefix

  Tr%conv_scale = 1.0
  if (present(convergence_scale)) then
    Tr%conv_scale = convergence_scale
  elseif (present(flux_scale)) then
    Tr%conv_scale = flux_scale
  endif

  Tr%diag_form = 1
  if (present(diag_form)) Tr%diag_form = diag_form

  Tr%t => tr_ptr

  if (present(registry_diags)) Tr%registry_diags = registry_diags

  if (present(ad_x)) then ; if (associated(ad_x)) Tr%ad_x => ad_x ; endif
  if (present(ad_y)) then ; if (associated(ad_y)) Tr%ad_y => ad_y ; endif
  if (present(df_x)) then ; if (associated(df_x)) Tr%df_x => df_x ; endif
  if (present(df_y)) then ; if (associated(df_y)) Tr%df_y => df_y ; endif
! if (present(OBC_inflow)) Tr%OBC_inflow_conc = OBC_inflow
! if (present(OBC_in_u)) then ; if (associated(OBC_in_u)) &
!                                   Tr%OBC_in_u => OBC_in_u ; endif
! if (present(OBC_in_v)) then ; if (associated(OBC_in_v)) &
!                                   Tr%OBC_in_v => OBC_in_v ; endif
  if (present(ad_2d_x)) then ; if (associated(ad_2d_x)) Tr%ad2d_x => ad_2d_x ; endif
  if (present(ad_2d_y)) then ; if (associated(ad_2d_y)) Tr%ad2d_y => ad_2d_y ; endif
  if (present(df_2d_x)) then ; if (associated(df_2d_x)) Tr%df2d_x => df_2d_x ; endif

  if (present(advection_xy)) then ; if (associated(advection_xy)) Tr%advection_xy => advection_xy ; endif

  if (present(restart_CS)) then ; if (associated(restart_CS)) then
    ! Register this tracer to be read from and written to restart files.
    mand = .true. ; if (present(mandatory)) mand = mandatory

    call register_restart_field(tr_ptr, Tr%name, mand, restart_CS, &
                                longname=Tr%longname, units=Tr%units)
  endif ; endif

end subroutine register_tracer


!> This subroutine locks the tracer registry to prevent the addition of more
!! tracers.  After locked=.true., can then register common diagnostics.
subroutine lock_tracer_registry(Reg)
  type(tracer_registry_type), pointer    :: Reg    !< pointer to the tracer registry

  if (.not. associated(Reg)) call MOM_error(WARNING, &
    "lock_tracer_registry called with an unassociated registry.")

  Reg%locked = .True.

end subroutine lock_tracer_registry

!> register_tracer_diagnostics does a set of register_diag_field calls for any previously
!! registered in a tracer registry with a value of registry_diags set to .true.
subroutine register_tracer_diagnostics(Reg, h, Time, diag, G, GV, US, use_ALE)
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in) :: h    !< Layer thicknesses
  type(time_type),            intent(in) :: Time !< current model time
  type(diag_ctrl),            intent(in) :: diag !< structure to regulate diagnostic output
  logical,                    intent(in) :: use_ALE !< If true active diagnostics that only
                                                 !! apply to ALE configurations

  character(len=24) :: name     ! A variable's name in a NetCDF file.
  character(len=24) :: shortnm  ! A shortened version of a variable's name for
                                ! creating additional diagnostics.
  character(len=72) :: longname ! The long name of that tracer variable.
  character(len=72) :: flux_longname ! The tracer name in the long names of fluxes.
  character(len=48) :: units    ! The dimensions of the tracer.
  character(len=48) :: flux_units ! The units for fluxes, either
                                ! [units] m3 s-1 or [units] kg s-1.
  character(len=48) :: conv_units ! The units for flux convergences, either
                                ! [units] m2 s-1 or [units] kg s-1.
  character(len=48) :: unit2    ! The dimensions of the tracer squared
  character(len=72)  :: cmorname ! The CMOR name of this tracer.
  character(len=120) :: cmor_longname ! The CMOR long name of that variable.
  character(len=120) :: var_lname      ! A temporary longname for a diagnostic.
  character(len=120) :: cmor_var_lname ! The temporary CMOR long name for a diagnostic
  character(len=72)  :: cmor_varname ! The temporary CMOR name for a diagnostic
  type(tracer_type), pointer :: Tr=>NULL()
  integer :: i, j, k, is, ie, js, je, nz, m, m2, nTr_in
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not. associated(Reg)) call MOM_error(FATAL, "register_tracer_diagnostics: "//&
       "register_tracer must be called before register_tracer_diagnostics")

  nTr_in = Reg%ntr

  do m=1,nTr_in ; if (Reg%Tr(m)%registry_diags) then
    Tr => Reg%Tr(m)
!    call query_vardesc(Tr%vd, name, units=units, longname=longname, &
!                       cmor_field_name=cmorname, cmor_longname=cmor_longname, &
!                       caller="register_tracer_diagnostics")
    name = Tr%name ; units=adjustl(Tr%units) ; longname = Tr%longname
    cmorname = Tr%cmor_name ; cmor_longname = Tr%cmor_longname
    shortnm = Tr%flux_nameroot
    flux_longname = Tr%flux_longname
    if (len_trim(cmor_longname) == 0) cmor_longname = longname

    if (len_trim(Tr%flux_units) > 0) then ; flux_units = Tr%flux_units
    elseif (GV%Boussinesq) then ; flux_units = trim(units)//" m3 s-1"
    else ; flux_units = trim(units)//" kg s-1" ; endif

    if (len_trim(Tr%conv_units) > 0) then ; conv_units = Tr%conv_units
    elseif (GV%Boussinesq) then ; conv_units = trim(units)//" m s-1"
    else ; conv_units = trim(units)//" kg m-2 s-1" ; endif

    if (len_trim(cmorname) == 0) then
      Tr%id_tr = register_diag_field("ocean_model", trim(name), diag%axesTL, &
        Time, trim(longname), trim(units))
    else
      Tr%id_tr = register_diag_field("ocean_model", trim(name), diag%axesTL, &
        Time, trim(longname), trim(units), cmor_field_name=cmorname, &
        cmor_long_name=cmor_longname, cmor_units=Tr%cmor_units, &
        cmor_standard_name=cmor_long_std(cmor_longname))
    endif
    if (Tr%diag_form == 1) then
      Tr%id_adx = register_diag_field("ocean_model", trim(shortnm)//"_adx", &
          diag%axesCuL, Time, trim(flux_longname)//" advective zonal flux" , &
          trim(flux_units), v_extensive = .true., y_cell_method = 'sum', &
          conversion=Tr%flux_scale*(US%L_to_m**2)*US%s_to_T)
      Tr%id_ady = register_diag_field("ocean_model", trim(shortnm)//"_ady", &
          diag%axesCvL, Time, trim(flux_longname)//" advective meridional flux" , &
          trim(flux_units), v_extensive = .true., x_cell_method = 'sum', &
          conversion=Tr%flux_scale*(US%L_to_m**2)*US%s_to_T)
      Tr%id_dfx = register_diag_field("ocean_model", trim(shortnm)//"_dfx", &
          diag%axesCuL, Time, trim(flux_longname)//" diffusive zonal flux" , &
          trim(flux_units), v_extensive = .true., y_cell_method = 'sum', &
          conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T)
      Tr%id_dfy = register_diag_field("ocean_model", trim(shortnm)//"_dfy", &
          diag%axesCvL, Time, trim(flux_longname)//" diffusive zonal flux" , &
          trim(flux_units), v_extensive = .true., x_cell_method = 'sum', &
          conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T)
      Tr%id_lbd_dfx = register_diag_field("ocean_model", trim(shortnm)//"_lbd_diffx", &
          diag%axesCuL, Time, trim(flux_longname)//" diffusive zonal flux from the lateral boundary diffusion "//&
          "scheme", trim(flux_units), v_extensive = .true., y_cell_method = 'sum', &
          conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T)
      Tr%id_lbd_dfy = register_diag_field("ocean_model", trim(shortnm)//"_lbd_diffy", &
          diag%axesCvL, Time, trim(flux_longname)//" diffusive meridional flux from the lateral boundary diffusion "//&
          "scheme", trim(flux_units), v_extensive = .true., x_cell_method = 'sum', &
          conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T)
    else
      Tr%id_adx = register_diag_field("ocean_model", trim(shortnm)//"_adx", &
          diag%axesCuL, Time, "Advective (by residual mean) Zonal Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=Tr%flux_scale*(US%L_to_m**2)*US%s_to_T, y_cell_method = 'sum')
      Tr%id_ady = register_diag_field("ocean_model", trim(shortnm)//"_ady", &
          diag%axesCvL, Time, "Advective (by residual mean) Meridional Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=Tr%flux_scale*(US%L_to_m**2)*US%s_to_T, x_cell_method = 'sum')
      Tr%id_dfx = register_diag_field("ocean_model", trim(shortnm)//"_diffx", &
          diag%axesCuL, Time, "Diffusive Zonal Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
          y_cell_method='sum')
      Tr%id_dfy = register_diag_field("ocean_model", trim(shortnm)//"_diffy", &
          diag%axesCvL, Time, "Diffusive Meridional Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
          x_cell_method='sum')
      Tr%id_lbd_dfx = register_diag_field("ocean_model", trim(shortnm)//"_lbd_diffx", &
          diag%axesCuL, Time, "Lateral Boundary Diffusive Zonal Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
          y_cell_method='sum')
      Tr%id_lbd_dfy = register_diag_field("ocean_model", trim(shortnm)//"_lbd_diffy", &
          diag%axesCvL, Time, "Lateral Boundary Diffusive Meridional Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
          x_cell_method='sum')
    endif
    if (Tr%id_adx > 0) call safe_alloc_ptr(Tr%ad_x,IsdB,IedB,jsd,jed,nz)
    if (Tr%id_ady > 0) call safe_alloc_ptr(Tr%ad_y,isd,ied,JsdB,JedB,nz)
    if (Tr%id_dfx > 0) call safe_alloc_ptr(Tr%df_x,IsdB,IedB,jsd,jed,nz)
    if (Tr%id_dfy > 0) call safe_alloc_ptr(Tr%df_y,isd,ied,JsdB,JedB,nz)
    if (Tr%id_lbd_dfx > 0) call safe_alloc_ptr(Tr%lbd_dfx,IsdB,IedB,jsd,jed,nz)
    if (Tr%id_lbd_dfy > 0) call safe_alloc_ptr(Tr%lbd_dfy,isd,ied,JsdB,JedB,nz)

    Tr%id_adx_2d = register_diag_field("ocean_model", trim(shortnm)//"_adx_2d", &
        diag%axesCu1, Time, &
        "Vertically Integrated Advective Zonal Flux of "//trim(flux_longname), &
        flux_units, conversion=Tr%flux_scale*(US%L_to_m**2)*US%s_to_T, y_cell_method = 'sum')
    Tr%id_ady_2d = register_diag_field("ocean_model", trim(shortnm)//"_ady_2d", &
        diag%axesCv1, Time, &
        "Vertically Integrated Advective Meridional Flux of "//trim(flux_longname), &
        flux_units, conversion=Tr%flux_scale*(US%L_to_m**2)*US%s_to_T, x_cell_method = 'sum')
    Tr%id_dfx_2d = register_diag_field("ocean_model", trim(shortnm)//"_diffx_2d", &
        diag%axesCu1, Time, &
        "Vertically Integrated Diffusive Zonal Flux of "//trim(flux_longname), &
        flux_units, conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
        y_cell_method='sum')
    Tr%id_dfy_2d = register_diag_field("ocean_model", trim(shortnm)//"_diffy_2d", &
        diag%axesCv1, Time, &
        "Vertically Integrated Diffusive Meridional Flux of "//trim(flux_longname), &
        flux_units, conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
        x_cell_method='sum')
    Tr%id_lbd_bulk_dfx = register_diag_field("ocean_model", trim(shortnm)//"_lbd_bulk_diffx", &
        diag%axesCu1, Time, &
        "Total Bulk Diffusive Zonal Flux of "//trim(flux_longname), &
        flux_units, conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
        y_cell_method='sum')
    Tr%id_lbd_bulk_dfy = register_diag_field("ocean_model", trim(shortnm)//"_lbd_bulk_diffy", &
        diag%axesCv1, Time, &
        "Total Bulk Diffusive Meridional Flux of "//trim(flux_longname), &
        flux_units, conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
        x_cell_method='sum')
    Tr%id_lbd_dfx_2d = register_diag_field("ocean_model", trim(shortnm)//"_lbd_diffx_2d", &
          diag%axesCu1, Time, "Vertically-integrated zonal diffusive flux from the lateral boundary diffusion "//&
          "scheme for "//trim(flux_longname), flux_units, conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
          y_cell_method = 'sum')
    Tr%id_lbd_dfy_2d = register_diag_field("ocean_model", trim(shortnm)//"_lbd_diffy_2d", &
          diag%axesCv1, Time, "Vertically-integrated meridional diffusive flux from the lateral boundary diffusion "//&
          "scheme for "//trim(flux_longname), flux_units, conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
           x_cell_method = 'sum')

    if (Tr%id_adx_2d > 0) call safe_alloc_ptr(Tr%ad2d_x,IsdB,IedB,jsd,jed)
    if (Tr%id_ady_2d > 0) call safe_alloc_ptr(Tr%ad2d_y,isd,ied,JsdB,JedB)
    if (Tr%id_dfx_2d > 0) call safe_alloc_ptr(Tr%df2d_x,IsdB,IedB,jsd,jed)
    if (Tr%id_dfy_2d > 0) call safe_alloc_ptr(Tr%df2d_y,isd,ied,JsdB,JedB)
    if (Tr%id_lbd_bulk_dfx > 0) call safe_alloc_ptr(Tr%lbd_bulk_df_x,IsdB,IedB,jsd,jed)
    if (Tr%id_lbd_bulk_dfy > 0) call safe_alloc_ptr(Tr%lbd_bulk_df_y,isd,ied,JsdB,JedB)
    if (Tr%id_lbd_dfx_2d > 0) call safe_alloc_ptr(Tr%lbd_dfx_2d,IsdB,IedB,jsd,jed)
    if (Tr%id_lbd_dfy_2d > 0) call safe_alloc_ptr(Tr%lbd_dfy_2d,isd,ied,JsdB,JedB)

    Tr%id_adv_xy = register_diag_field('ocean_model', trim(shortnm)//"_advection_xy", &
        diag%axesTL, Time, &
        'Horizontal convergence of residual mean advective fluxes of '//&
        trim(lowercase(flux_longname)), conv_units, v_extensive=.true., &
        conversion=Tr%conv_scale*US%s_to_T)
    Tr%id_adv_xy_2d = register_diag_field('ocean_model', trim(shortnm)//"_advection_xy_2d", &
        diag%axesT1, Time, &
        'Vertical sum of horizontal convergence of residual mean advective fluxes of '//&
        trim(lowercase(flux_longname)), conv_units, conversion=Tr%conv_scale*US%s_to_T)
    if ((Tr%id_adv_xy > 0) .or. (Tr%id_adv_xy_2d > 0)) &
      call safe_alloc_ptr(Tr%advection_xy,isd,ied,jsd,jed,nz)

    Tr%id_tendency = register_diag_field('ocean_model', trim(shortnm)//'_tendency', &
        diag%axesTL, Time, &
        'Net time tendency for '//trim(lowercase(longname)), trim(units)//' s-1', conversion=US%s_to_T)

    if (Tr%id_tendency > 0) then
      call safe_alloc_ptr(Tr%t_prev,isd,ied,jsd,jed,nz)
      do k=1,nz ; do j=js,je ; do i=is,ie
        Tr%t_prev(i,j,k) = Tr%t(i,j,k)
      enddo ; enddo ; enddo
    endif

    ! Neutral/Lateral diffusion convergence tendencies
    if (Tr%diag_form == 1) then
      Tr%id_dfxy_cont = register_diag_field("ocean_model", trim(shortnm)//'_dfxy_cont_tendency', &
          diag%axesTL, Time, "Neutral diffusion tracer content tendency for "//trim(shortnm), &
          conv_units, conversion=Tr%conv_scale*US%s_to_T, x_cell_method='sum', y_cell_method='sum', v_extensive=.true.)

      Tr%id_dfxy_cont_2d = register_diag_field("ocean_model", trim(shortnm)//'_dfxy_cont_tendency_2d', &
          diag%axesT1, Time, "Depth integrated neutral diffusion tracer content "//&
          "tendency for "//trim(shortnm), conv_units, conversion=Tr%conv_scale*US%s_to_T, &
          x_cell_method='sum', y_cell_method= 'sum')

      Tr%id_lbdxy_cont = register_diag_field("ocean_model", trim(shortnm)//'_lbdxy_cont_tendency', &
          diag%axesTL, Time, "Lateral diffusion tracer content tendency for "//trim(shortnm), &
          conv_units, conversion=Tr%conv_scale*US%s_to_T, x_cell_method='sum', y_cell_method='sum', v_extensive=.true.)

      Tr%id_lbdxy_cont_2d = register_diag_field("ocean_model", trim(shortnm)//'_lbdxy_cont_tendency_2d', &
          diag%axesT1, Time, "Depth integrated lateral diffusion tracer content "//&
          "tendency for "//trim(shortnm), conv_units, conversion=Tr%conv_scale*US%s_to_T, &
          x_cell_method='sum', y_cell_method= 'sum')
    else
      cmor_var_lname = 'Tendency of '//trim(lowercase(cmor_longname))//&
           ' expressed as '//trim(lowercase(flux_longname))//&
           ' content due to parameterized mesoscale neutral diffusion'
      Tr%id_dfxy_cont = register_diag_field("ocean_model", trim(shortnm)//'_dfxy_cont_tendency', &
          diag%axesTL, Time, "Neutral diffusion tracer content tendency for "//trim(shortnm), &
          conv_units, conversion=Tr%conv_scale*US%s_to_T, cmor_field_name = trim(Tr%cmor_tendprefix)//'pmdiff', &
          cmor_long_name = trim(cmor_var_lname), cmor_standard_name = trim(cmor_long_std(cmor_var_lname)), &
          x_cell_method = 'sum', y_cell_method = 'sum', v_extensive = .true.)

      cmor_var_lname = 'Tendency of '//trim(lowercase(cmor_longname))//' expressed as '//&
                       trim(lowercase(flux_longname))//' content due to parameterized mesoscale neutral diffusion'
      Tr%id_dfxy_cont_2d = register_diag_field("ocean_model", trim(shortnm)//'_dfxy_cont_tendency_2d', &
          diag%axesT1, Time, "Depth integrated neutral diffusion tracer "//&
          "content tendency for "//trim(shortnm), conv_units, &
          conversion=Tr%conv_scale*US%s_to_T, cmor_field_name=trim(Tr%cmor_tendprefix)//'pmdiff_2d', &
          cmor_long_name=trim(cmor_var_lname), cmor_standard_name=trim(cmor_long_std(cmor_var_lname)), &
          x_cell_method='sum', y_cell_method='sum')

      Tr%id_lbdxy_cont = register_diag_field("ocean_model", trim(shortnm)//'_lbdxy_cont_tendency', &
          diag%axesTL, Time, "Lateral diffusion tracer content tendency for "//trim(shortnm), &
          conv_units, conversion=Tr%conv_scale*US%s_to_T, &
          x_cell_method = 'sum', y_cell_method = 'sum', v_extensive = .true.)

      Tr%id_lbdxy_cont_2d = register_diag_field("ocean_model", trim(shortnm)//'_lbdxy_cont_tendency_2d', &
          diag%axesT1, Time, "Depth integrated lateral diffusion tracer "//&
          "content tendency for "//trim(shortnm), conv_units, &
          conversion=Tr%conv_scale*US%s_to_T, x_cell_method='sum', y_cell_method='sum')
    endif
    Tr%id_dfxy_conc = register_diag_field("ocean_model", trim(shortnm)//'_dfxy_conc_tendency', &
        diag%axesTL, Time, "Neutral diffusion tracer concentration tendency for "//trim(shortnm), &
        trim(units)//' s-1', conversion=US%s_to_T)

    Tr%id_lbdxy_conc = register_diag_field("ocean_model", trim(shortnm)//'_lbdxy_conc_tendency', &
        diag%axesTL, Time, "Lateral diffusion tracer concentration tendency for "//trim(shortnm), &
        trim(units)//' s-1', conversion=US%s_to_T)

    var_lname = "Net time tendency for "//lowercase(flux_longname)
    if (len_trim(Tr%cmor_tendprefix) == 0) then
      Tr%id_trxh_tendency = register_diag_field('ocean_model', trim(shortnm)//'h_tendency', &
          diag%axesTL, Time, var_lname, conv_units, v_extensive=.true., &
          conversion=Tr%conv_scale*US%s_to_T)
      Tr%id_trxh_tendency_2d = register_diag_field('ocean_model', trim(shortnm)//'h_tendency_2d', &
          diag%axesT1, Time, "Vertical sum of "//trim(lowercase(var_lname)), conv_units, &
          conversion=Tr%conv_scale*US%s_to_T)
    else
      cmor_var_lname = "Tendency of "//trim(cmor_longname)//" Expressed as "//&
                        trim(flux_longname)//" Content"
      Tr%id_trxh_tendency = register_diag_field('ocean_model', trim(shortnm)//'h_tendency', &
          diag%axesTL, Time, var_lname, conv_units, &
          cmor_field_name=trim(Tr%cmor_tendprefix)//"tend", &
          cmor_standard_name=cmor_long_std(cmor_var_lname), cmor_long_name=cmor_var_lname, &
          v_extensive=.true., conversion=Tr%conv_scale*US%s_to_T)
      cmor_var_lname = trim(cmor_var_lname)//" Vertical Sum"
      Tr%id_trxh_tendency_2d = register_diag_field('ocean_model', trim(shortnm)//'h_tendency_2d', &
          diag%axesT1, Time, "Vertical sum of "//trim(lowercase(var_lname)), conv_units, &
          cmor_field_name=trim(Tr%cmor_tendprefix)//"tend_2d", &
          cmor_standard_name=cmor_long_std(cmor_var_lname), cmor_long_name=cmor_var_lname, &
          conversion=Tr%conv_scale*US%s_to_T)
    endif
    if ((Tr%id_trxh_tendency > 0) .or. (Tr%id_trxh_tendency_2d > 0)) then
      call safe_alloc_ptr(Tr%Trxh_prev,isd,ied,jsd,jed,nz)
      do k=1,nz ; do j=js,je ; do i=is,ie
        Tr%Trxh_prev(i,j,k) = Tr%t(i,j,k) * h(i,j,k)
      enddo ; enddo ; enddo
    endif

    ! Vertical regridding/remapping tendencies
    if (use_ALE .and. Tr%remap_tr) then
      var_lname = "Vertical remapping tracer concentration tendency for "//trim(Reg%Tr(m)%name)
      Tr%id_remap_conc= register_diag_field('ocean_model',                          &
        trim(Tr%flux_nameroot)//'_tendency_vert_remap', diag%axesTL, Time, var_lname, &
        trim(units)//' s-1', conversion=US%s_to_T)

      var_lname = "Vertical remapping tracer content tendency for "//trim(Reg%Tr(m)%flux_longname)
      Tr%id_remap_cont = register_diag_field('ocean_model', &
        trim(Tr%flux_nameroot)//'h_tendency_vert_remap',         &
        diag%axesTL, Time, var_lname, flux_units, v_extensive=.true., conversion=Tr%conv_scale*US%s_to_T)

      var_lname = "Vertical sum of vertical remapping tracer content tendency for "//&
                  trim(Reg%Tr(m)%flux_longname)
      Tr%id_remap_cont_2d = register_diag_field('ocean_model', &
        trim(Tr%flux_nameroot)//'h_tendency_vert_remap_2d',         &
        diag%axesT1, Time, var_lname, flux_units, conversion=Tr%conv_scale*US%s_to_T)

    endif

    if (use_ALE .and. (Reg%ntr<MAX_FIELDS_) .and. Tr%remap_tr) then
      unit2 = trim(units)//"2"
      if (index(units(1:len_trim(units))," ") > 0) unit2 = "("//trim(units)//")2"
      Tr%id_tr_vardec = register_diag_field('ocean_model', trim(shortnm)//"_vardec", diag%axesTL, &
        Time, "ALE variance decay for "//lowercase(longname), trim(unit2)//" s-1", conversion=US%s_to_T)
      if (Tr%id_tr_vardec > 0) then
        ! Set up a new tracer for this tracer squared
        m2 = Reg%ntr+1
        Tr%ind_tr_squared = m2
        call safe_alloc_ptr(Reg%Tr(m2)%t,isd,ied,jsd,jed,nz) ; Reg%Tr(m2)%t(:,:,:) = 0.0
        Reg%Tr(m2)%name = trim(shortnm)//"2"
        Reg%Tr(m2)%longname = "Squared "//trim(longname)
        Reg%Tr(m2)%units = unit2
        Reg%Tr(m2)%registry_diags = .false.
        Reg%Tr(m2)%ind_tr_squared = -1
        ! Augment the total number of tracers, including the squared tracers.
        Reg%ntr = Reg%ntr + 1
      endif
    endif

  endif ; enddo

end subroutine register_tracer_diagnostics

subroutine preALE_tracer_diagnostics(Reg, G, GV)
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< ocean vertical grid structure

  integer :: i, j, k, is, ie, js, je, nz, m, m2
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  do m=1,Reg%ntr ; if (Reg%Tr(m)%ind_tr_squared > 0) then
    m2 = Reg%Tr(m)%ind_tr_squared
  ! Update squared quantities
    do k=1,nz ; do j=js,je ; do i=is,ie
      Reg%Tr(m2)%T(i,j,k) = Reg%Tr(m)%T(i,j,k)**2
    enddo ; enddo ; enddo
  endif ; enddo

end subroutine preALE_tracer_diagnostics

subroutine postALE_tracer_diagnostics(Reg, G, GV, diag, dt)
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< ocean vertical grid structure
  type(diag_ctrl),            intent(in) :: diag !< regulates diagnostic output
  real,                       intent(in) :: dt   !< total time interval for these diagnostics [T ~> s]

  real    :: work(SZI_(G),SZJ_(G),SZK_(G))
  real    :: Idt ! The inverse of the time step [T-1 ~> s-1]
  integer :: i, j, k, is, ie, js, je, nz, m, m2
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! The "if" is to avoid NaNs if the diagnostic is called for a zero length interval
  Idt = 0.0 ; if (dt /= 0.0) Idt = 1.0 / dt

  do m=1,Reg%ntr ; if (Reg%Tr(m)%id_tr_vardec > 0) then
    m2 = Reg%Tr(m)%ind_tr_squared
    if (m2 < 1) call MOM_error(FATAL, "Bad value of Tr%ind_tr_squared for "//trim(Reg%Tr(m)%name))
  ! Update squared quantities
    do k=1,nz ; do j=js,je ; do i=is,ie
      work(i,j,k) = (Reg%Tr(m2)%T(i,j,k) - Reg%Tr(m)%T(i,j,k)**2) * Idt
    enddo ; enddo ; enddo
    call post_data(Reg%Tr(m)%id_tr_vardec, work, diag)
  endif ; enddo

end subroutine postALE_tracer_diagnostics

!> post_tracer_diagnostics does post_data calls for any diagnostics that are
!! being handled via the tracer registry.
subroutine post_tracer_diagnostics(Reg, h, diag_prev, diag, G, GV, dt)
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in) :: h    !< Layer thicknesses
  type(diag_grid_storage),    intent(in) :: diag_prev !< Contains diagnostic grids from previous timestep
  type(diag_ctrl),            intent(inout) :: diag !< structure to regulate diagnostic output
  real,                       intent(in) :: dt   !< total time step for tracer updates [T ~> s]

  real    :: work3d(SZI_(G),SZJ_(G),SZK_(GV))
  real    :: work2d(SZI_(G),SZJ_(G))
  real    :: Idt ! The inverse of the time step [T-1 ~> s-1]
  type(tracer_type), pointer :: Tr=>NULL()
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  Idt = 0.; if (dt/=0.) Idt = 1.0 / dt ! The "if" is in case the diagnostic is called for a zero length interval

  ! Tendency diagnostics need to be posted on the grid from the last call to this routine
  call diag_save_grids(diag)
  call diag_copy_storage_to_diag(diag, diag_prev)
  do m=1,Reg%ntr ; if (Reg%Tr(m)%registry_diags) then
    Tr => Reg%Tr(m)
    if (Tr%id_tendency > 0) then
      work3d(:,:,:) = 0.0
      do k=1,nz ; do j=js,je ; do i=is,ie
        work3d(i,j,k)    = (Tr%t(i,j,k) - Tr%t_prev(i,j,k))*Idt
        tr%t_prev(i,j,k) =  Tr%t(i,j,k)
      enddo ; enddo ; enddo
      call post_data(Tr%id_tendency, work3d, diag, alt_h=diag_prev%h_state)
    endif
    if ((Tr%id_trxh_tendency > 0) .or. (Tr%id_trxh_tendency_2d > 0)) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        work3d(i,j,k)     = (Tr%t(i,j,k)*h(i,j,k) - Tr%Trxh_prev(i,j,k)) * Idt
        Tr%Trxh_prev(i,j,k) =  Tr%t(i,j,k) * h(i,j,k)
      enddo ; enddo ; enddo
      if (Tr%id_trxh_tendency > 0) call post_data(Tr%id_trxh_tendency, work3d, diag, alt_h=diag_prev%h_state)
      if (Tr%id_trxh_tendency_2d > 0) then
        work2d(:,:) = 0.0
        do k=1,nz ; do j=js,je ; do i=is,ie
          work2d(i,j) = work2d(i,j) + work3d(i,j,k)
        enddo ; enddo ; enddo
        call post_data(Tr%id_trxh_tendency_2d, work2d, diag)
      endif
    endif
  endif ; enddo
  call diag_restore_grids(diag)

end subroutine post_tracer_diagnostics

!> Post the advective and diffusive tendencies
subroutine post_tracer_transport_diagnostics(G, GV, Reg, h_diag, diag)
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in) :: h_diag !< Layer thicknesses on which to post fields
  type(diag_ctrl),            intent(in) :: diag !< structure to regulate diagnostic output

  integer :: i, j, k, is, ie, js, je, nz, m
  real    :: work2d(SZI_(G),SZJ_(G))
  type(tracer_type), pointer :: Tr=>NULL()

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  do m=1,Reg%ntr ; if (Reg%Tr(m)%registry_diags) then
    Tr => Reg%Tr(m)
    if (Tr%id_tr > 0) call post_data(Tr%id_tr, Tr%t, diag)
    if (Tr%id_adx > 0) call post_data(Tr%id_adx, Tr%ad_x, diag, alt_h=h_diag)
    if (Tr%id_ady > 0) call post_data(Tr%id_ady, Tr%ad_y, diag, alt_h=h_diag)
    if (Tr%id_dfx > 0) call post_data(Tr%id_dfx, Tr%df_x, diag, alt_h=h_diag)
    if (Tr%id_dfy > 0) call post_data(Tr%id_dfy, Tr%df_y, diag, alt_h=h_diag)
    if (Tr%id_adx_2d > 0) call post_data(Tr%id_adx_2d, Tr%ad2d_x, diag)
    if (Tr%id_ady_2d > 0) call post_data(Tr%id_ady_2d, Tr%ad2d_y, diag)
    if (Tr%id_dfx_2d > 0) call post_data(Tr%id_dfx_2d, Tr%df2d_x, diag)
    if (Tr%id_dfy_2d > 0) call post_data(Tr%id_dfy_2d, Tr%df2d_y, diag)
    if (Tr%id_adv_xy > 0) call post_data(Tr%id_adv_xy, Tr%advection_xy, diag, alt_h=h_diag)
    if (Tr%id_adv_xy_2d > 0) then
      work2d(:,:) = 0.0
      do k=1,nz ; do j=js,je ; do i=is,ie
        work2d(i,j) = work2d(i,j) + Tr%advection_xy(i,j,k)
      enddo ; enddo ; enddo
      call post_data(Tr%id_adv_xy_2d, work2d, diag)
    endif
  endif ; enddo

end subroutine post_tracer_transport_diagnostics

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
  character(len=*),                         intent(in) :: mesg !< message that appears on the chksum lines
  type(ocean_grid_type),                    intent(in) :: G    !< ocean grid structure
  type(tracer_type), dimension(:),          intent(in) :: Tr   !< array of all of registered tracers
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h    !< Layer thicknesses
  integer,                                  intent(in) :: ntr  !< number of registered tracers

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: tr_inv !< Tracer inventory
  real :: total_inv
  integer :: is, ie, js, je, nz
  integer :: i, j, k, m

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  do m=1,ntr
    do k=1,nz ; do j=js,je ; do i=is,ie
      tr_inv(i,j,k) = Tr(m)%t(i,j,k)*h(i,j,k)*G%US%L_to_m**2*G%areaT(i,j)*G%mask2dT(i,j)
    enddo ; enddo ; enddo
    total_inv = reproducing_sum(tr_inv, is+(1-G%isd), ie+(1-G%isd), js+(1-G%jsd), je+(1-G%jsd))
    if (is_root_pe()) write(0,'(A,1X,A5,1X,ES25.16,1X,A)') "h-point: inventory", Tr(m)%name, total_inv, mesg
  enddo

end subroutine MOM_tracer_chkinv

!> Find a tracer in the tracer registry by name.
subroutine tracer_name_lookup(Reg, tr_ptr, name)
  type(tracer_registry_type), pointer    :: Reg     !< pointer to tracer registry
  type(tracer_type), pointer             :: tr_ptr  !< target or pointer to the tracer array
  character(len=32), intent(in)          :: name    !< tracer name

  integer n
  do n=1,Reg%ntr
    if (lowercase(Reg%Tr(n)%name) == lowercase(name)) tr_ptr => Reg%Tr(n)
  enddo

end subroutine tracer_name_lookup

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
  type(tracer_registry_type), pointer :: Reg !< The tracer registry that will be deallocated
  if (associated(Reg)) deallocate(Reg)
end subroutine tracer_registry_end

end module MOM_tracer_registry
