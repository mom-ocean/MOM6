!> This module contains the tracer_type and tracer_registry_type
module MOM_tracer_types

implicit none ; private

#include <MOM_memory.h>

!> The tracer type
type, public :: tracer_type

  real, dimension(:,:,:), pointer :: t              => NULL() !< tracer concentration array [CU ~> conc]
! real                            :: OBC_inflow_conc=  0.0    !< tracer concentration for generic inflows [CU ~> conc]
! real, dimension(:,:,:), pointer :: OBC_in_u       => NULL() !< structured values for flow into the domain
!                                                             !! specified in OBCs through u-face of cell
! real, dimension(:,:,:), pointer :: OBC_in_v       => NULL() !< structured values for flow into the domain
!                                                             !! specified in OBCs through v-face of cell

  real, dimension(:,:,:), pointer :: ad_x           => NULL() !< diagnostic array for x-advective tracer flux
                                                              !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:,:), pointer :: ad_y           => NULL() !< diagnostic array for y-advective tracer flux
                                                              !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   pointer :: ad2d_x         => NULL() !< diagnostic vertical sum x-advective tracer flux
                                                              !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   pointer :: ad2d_y         => NULL() !< diagnostic vertical sum y-advective tracer flux
                                                              !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]

  real, dimension(:,:,:), pointer :: df_x           => NULL() !< diagnostic array for x-diffusive tracer flux
                                                              !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:,:), pointer :: df_y           => NULL() !< diagnostic array for y-diffusive tracer flux
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:,:), pointer :: hbd_dfx       => NULL()  !< diagnostic array for x-diffusive tracer flux
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:,:), pointer :: hbd_dfy       => NULL()  !< diagnostic array for y-diffusive tracer flux
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   pointer :: hbd_dfx_2d    => NULL()  !< diagnostic array for x-diffusive tracer flux
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   pointer :: hbd_dfy_2d    => NULL()  !< diagnostic array for y-diffusive tracer flux
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
                                                              !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   pointer :: df2d_x         => NULL() !< diagnostic vertical sum x-diffusive flux
                                                              !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   pointer :: df2d_y         => NULL() !< diagnostic vertical sum y-diffusive flux
                                                              !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
!  real, dimension(:,:),   pointer :: df2d_conc_x    => NULL() !< diagnostic vertical sum x-diffusive content flux
!                                                              !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
!  real, dimension(:,:),   pointer :: df2d_conc_y    => NULL() !< diagnostic vertical sum y-diffusive content flux
!                                                              !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]

  real, dimension(:,:,:), pointer :: advection_xy   => NULL() !< convergence of lateral advective tracer fluxes
                                                              !! [CU H T-1 ~> conc m s-1 or conc kg m-2 s-1]
!  real, dimension(:,:,:), pointer :: diff_cont_xy   => NULL() !< convergence of lateral diffusive tracer fluxes
!                                                              !! [CU H T-1 ~> conc m s-1 or conc kg m-2 s-1]
!  real, dimension(:,:,:), pointer :: diff_conc_xy   => NULL() !< convergence of lateral diffusive tracer fluxes
!                                                              !! expressed as a change in concentration
!                                                              !! [CU T-1 ~> conc s-1]
  real, dimension(:,:,:), pointer :: t_prev         => NULL() !< tracer concentration array at a previous
                                                              !! timestep used for diagnostics [CU ~> conc]
  real, dimension(:,:,:), pointer :: Trxh_prev      => NULL() !< layer integrated tracer concentration array
                                                              !! at a previous timestep used for diagnostics
                                                              !! [CU H ~> conc m or conc kg m-2]

  character(len=32)               :: name                     !< tracer name used for diagnostics and error messages
  character(len=64)               :: units                    !< Physical dimensions of the tracer concentration
  character(len=240)              :: longname                 !< Long name of the variable
!  type(vardesc), pointer          :: vd             => NULL() !< metadata describing the tracer
  logical                         :: registry_diags = .false. !< If true, use the registry to set up the
                                                              !! diagnostics associated with this tracer.
  real                            :: conc_underflow = 0.0     !< A magnitude of tracer concentrations below
                                                              !! which values should be set to 0. [CU ~> conc]
  real                            :: conc_scale = 1.0         !< A scaling factor used to convert the concentrations
                                                              !! of this tracer to its desired units [CU conc-1 ~> 1]
  character(len=64)               :: cmor_name                !< CMOR name of this tracer
  character(len=64)               :: cmor_units               !< CMOR physical dimensions of the tracer
  character(len=240)              :: cmor_longname            !< CMOR long name of the tracer
  character(len=32)               :: flux_nameroot = ""       !< Short tracer name snippet used construct the
                                                              !! names of flux diagnostics.
  character(len=64)               :: flux_longname = ""       !< A word or phrase used construct the long
                                                              !! names of flux diagnostics.
  real                            :: flux_scale = 1.0         !< A scaling factor used to convert the fluxes
                                                              !! of this tracer to its desired units,
                                                              !! including a factor compensating for H scaling.
                                                              !! [conc m CU-1 H-1 ~> 1] or [conc kg m-2 CU-1 H-1 ~> 1]
  character(len=48)               :: flux_units = ""          !< The units for fluxes of this variable.
  character(len=48)               :: conv_units = ""          !< The units for the flux convergence of this tracer.
  real                            :: conv_scale = 1.0         !< A scaling factor used to convert the flux
                                                              !! convergence of this tracer to its desired units,
                                                              !! including a factor compensating for H scaling.
                                                              !! [conc m CU-1 H-1 ~> 1] or [conc kg m-2 CU-1 H-1 ~> 1]
  character(len=48)               :: cmor_tendprefix = ""     !< The CMOR variable prefix for tendencies of this
                                                              !! tracer, required because CMOR does not follow any
                                                              !! discernable pattern for these names.
  character(len=48)               :: net_surfflux_name = ""   !< Name to use for net_surfflux KPP diagnostic
  character(len=48)               :: NLT_budget_name = ""     !< Name to use for NLT_budget KPP diagnostic
  character(len=128)              :: net_surfflux_longname = ""   !< Long name to use for net_surfflux KPP diagnostic
  integer :: ind_tr_squared = -1 !< The tracer registry index for the square of this tracer

  !### THESE CAPABILITIES HAVE NOT YET BEEN IMPLEMENTED.
  ! logical :: advect_tr = .true.       !< If true, this tracer should be advected
  ! logical :: hordiff_tr = .true.      !< If true, this tracer should experience epineutral diffusion
  ! logical :: kpp_nonlocal_tr = .true. !< if true, apply KPP nonlocal transport to this tracer before diffusion
  logical :: remap_tr = .true.        !< If true, this tracer should be vertically remapped
  integer :: advect_scheme = -1  !< flag for advection scheme

  integer :: diag_form = 1  !< An integer indicating which template is to be used to label diagnostics.
  !>@{ Diagnostic IDs
  integer :: id_tr = -1, id_tr_post_horzn = -1
  integer :: id_adx = -1, id_ady = -1, id_dfx = -1, id_dfy = -1
  integer :: id_hbd_dfx = -1, id_hbd_dfy = -1
  integer :: id_hbd_dfx_2d = -1, id_hbd_dfy_2d = -1
  integer :: id_adx_2d = -1, id_ady_2d = -1, id_dfx_2d = -1, id_dfy_2d = -1
  integer :: id_adv_xy = -1, id_adv_xy_2d = -1
  integer :: id_dfxy_cont = -1, id_dfxy_cont_2d = -1, id_dfxy_conc = -1
  integer :: id_hbdxy_cont = -1, id_hbdxy_cont_2d = -1, id_hbdxy_conc = -1
  integer :: id_remap_conc = -1, id_remap_cont = -1, id_remap_cont_2d = -1
  integer :: id_tendency = -1, id_trxh_tendency = -1, id_trxh_tendency_2d = -1
  integer :: id_tr_vardec = -1
  integer :: id_zint = -1, id_zint_100m = -1, id_surf = -1
  integer :: id_net_surfflux = -1, id_NLT_tendency = -1, id_NLT_budget = -1
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


end module MOM_tracer_types
