!> Calculate vertical diffusivity from all mixing processes
module MOM_set_diffusivity

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_bkgnd_mixing,        only : calculate_bkgnd_mixing, bkgnd_mixing_init, bkgnd_mixing_cs
use MOM_bkgnd_mixing,        only : bkgnd_mixing_end
use MOM_cpu_clock,           only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,           only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_CVMix_ddiff,         only : CVMix_ddiff_init, CVMix_ddiff_end, CVMix_ddiff_cs
use MOM_CVMix_ddiff,         only : compute_ddiff_coeffs
use MOM_CVMix_shear,         only : calculate_CVMix_shear, CVMix_shear_init, CVMix_shear_cs
use MOM_CVMix_shear,         only : CVMix_shear_end
use MOM_diag_mediator,       only : diag_ctrl, time_type
use MOM_diag_mediator,       only : post_data, register_diag_field
use MOM_debugging,           only : hchksum, uvchksum, Bchksum, hchksum_pair
use MOM_EOS,                 only : calculate_density, calculate_density_derivs, EOS_domain
use MOM_error_handler,       only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_error_handler,       only : callTree_showQuery
use MOM_error_handler,       only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,         only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,        only : forcing, optics_type
use MOM_full_convection,     only : full_convection
use MOM_grid,                only : ocean_grid_type
use MOM_internal_tides,      only : int_tide_CS, get_lowmode_loss
use MOM_intrinsic_functions, only : invcosh
use MOM_io,                  only : slasher, MOM_read_data
use MOM_isopycnal_slopes,    only : vert_fill_TS
use MOM_kappa_shear,         only : calculate_kappa_shear, kappa_shear_init, Kappa_shear_CS
use MOM_kappa_shear,         only : calc_kappa_shear_vertex, kappa_shear_at_vertex
use MOM_open_boundary,       only : ocean_OBC_type, OBC_segment_type, OBC_NONE
use MOM_open_boundary,       only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_string_functions,    only : uppercase
use MOM_tidal_mixing,        only : tidal_mixing_CS, calculate_tidal_mixing, tidal_mixing_h_amp
use MOM_tidal_mixing,        only : setup_tidal_diagnostics, post_tidal_diagnostics
use MOM_tidal_mixing,        only : tidal_mixing_init, tidal_mixing_end
use MOM_unit_scaling,        only : unit_scale_type
use MOM_variables,           only : thermo_var_ptrs, vertvisc_type, p3d
use MOM_verticalGrid,        only : verticalGrid_type
use user_change_diffusivity, only : user_change_diff, user_change_diff_init
use user_change_diffusivity, only : user_change_diff_end, user_change_diff_CS

implicit none ; private

#include <MOM_memory.h>

public set_diffusivity
public set_BBL_TKE
public set_diffusivity_init
public set_diffusivity_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> This control structure contains parameters for MOM_set_diffusivity.
type, public :: set_diffusivity_CS ; private
  logical :: debug           !< If true, write verbose checksums for debugging.

  logical :: bulkmixedlayer  !< If true, a refined bulk mixed layer is used with
                             !! GV%nk_rho_varies variable density mixed & buffer layers.
  real    :: FluxRi_max      !< The flux Richardson number where the stratification is
                             !! large enough that N2 > omega2 [nondim].  The full expression
                             !! for the Flux Richardson number is usually
                             !! FLUX_RI_MAX*N2/(N2+OMEGA2). The default is 0.2.
  logical :: bottomdraglaw   !< If true, the  bottom stress is calculated with a
                             !! drag law c_drag*|u|*u.
  logical :: BBL_mixing_as_max !<  If true, take the maximum of the diffusivity
                             !! from the BBL mixing and the other diffusivities.
                             !! Otherwise, diffusivities from the BBL_mixing is
                             !! added.
  logical :: use_LOTW_BBL_diffusivity !< If true, use simpler/less precise, BBL diffusivity.
  logical :: LOTW_BBL_use_omega !< If true, use simpler/less precise, BBL diffusivity.
  real    :: BBL_effic       !< efficiency with which the energy extracted
                             !! by bottom drag drives BBL diffusion [nondim]
  real    :: cdrag           !< quadratic drag coefficient [nondim]
  real    :: IMax_decay      !< inverse of a maximum decay scale for
                             !! bottom-drag driven turbulence [Z-1 ~> m-1].
  real    :: Kv              !< The interior vertical viscosity [Z2 T-1 ~> m2 s-1].
  real    :: Kd              !< interior diapycnal diffusivity [Z2 T-1 ~> m2 s-1].
  real    :: Kd_min          !< minimum diapycnal diffusivity [Z2 T-1 ~> m2 s-1].
  real    :: Kd_max          !< maximum increment for diapycnal diffusivity [Z2 T-1 ~> m2 s-1].
                             !! Set to a negative value to have no limit.
  real    :: Kd_add          !< uniform diffusivity added everywhere without
                             !! filtering or scaling [Z2 T-1 ~> m2 s-1].
  real    :: Kd_smooth       !< Vertical diffusivity used to interpolate more
                             !! sensible values of T & S into thin layers [Z2 T-1 ~> m2 s-1].
  type(diag_ctrl), pointer :: diag => NULL() !< structure to regulate diagnostic output timing

  logical :: limit_dissipation !< If enabled, dissipation is limited to be larger
                               !! than the following:
  real :: dissip_min    !< Minimum dissipation [R Z2 T-3 ~> W m-3]
  real :: dissip_N0     !< Coefficient a in minimum dissipation = a+b*N [R Z2 T-3 ~> W m-3]
  real :: dissip_N1     !< Coefficient b in minimum dissipation = a+b*N [R Z2 T-2 ~> J m-3]
  real :: dissip_N2     !< Coefficient c in minimum dissipation = c*N2 [R Z2 T-1 ~> J s m-3]
  real :: dissip_Kd_min !< Minimum Kd [Z2 T-1 ~> m2 s-1], with dissipation Rho0*Kd_min*N^2

  real :: omega         !< Earth's rotation frequency [T-1 ~> s-1]
  logical :: ML_radiation !< allow a fraction of TKE available from wind work
                          !! to penetrate below mixed layer base with a vertical
                          !! decay scale determined by the minimum of
                          !! (1) The depth of the mixed layer, or
                          !! (2) An Ekman length scale.
                          !! Energy available to drive mixing below the mixed layer is
                          !! given by E = ML_RAD_COEFF*MSTAR*USTAR**3.  Optionally, if
                          !! ML_rad_TKE_decay is true, this is further reduced by a factor
                          !! of exp(-h_ML*Idecay_len_TkE), where Idecay_len_TKE is
                          !! calculated the same way as in the mixed layer code.
                          !! The diapycnal diffusivity is KD(k) = E/(N2(k)+OMEGA2),
                          !! where N2 is the squared buoyancy frequency [T-2 ~> s-2] and OMEGA2
                          !! is the rotation rate of the earth squared.
  real :: ML_rad_kd_max   !< Maximum diapycnal diffusivity due to turbulence
                          !! radiated from the base of the mixed layer [Z2 T-1 ~> m2 s-1].
  real :: ML_rad_efold_coeff  !< non-dim coefficient to scale penetration depth
  real :: ML_rad_coeff        !< coefficient, which scales MSTAR*USTAR^3 to
                              !! obtain energy available for mixing below
                              !! mixed layer base [nondim]
  logical :: ML_rad_bug       !< If true use code with a bug that reduces the energy available
                              !! in the transition layer by a factor of the inverse of the energy
                              !! deposition lenthscale (in m).
  logical :: ML_rad_TKE_decay !< If true, apply same exponential decay
                              !! to ML_rad as applied to the other surface
                              !! sources of TKE in the mixed layer code.
  real    :: ustar_min        !< A minimum value of ustar to avoid numerical
                              !! problems [Z T-1 ~> m s-1].  If the value is small enough,
                              !! this parameter should not affect the solution.
  real    :: TKE_decay        !< ratio of natural Ekman depth to TKE decay scale [nondim]
  real    :: mstar            !< ratio of friction velocity cubed to
                              !! TKE input to the mixed layer [nondim]
  logical :: ML_use_omega     !< If true, use absolute rotation rate instead
                              !! of the vertical component of rotation when
                              !! setting the decay scale for mixed layer turbulence.
  real    :: ML_omega_frac    !<   When setting the decay scale for turbulence, use
                              !! this fraction of the absolute rotation rate blended
                              !! with the local value of f, as f^2 ~= (1-of)*f^2 + of*4*omega^2.
  logical :: user_change_diff !< If true, call user-defined code to change diffusivity.
  logical :: useKappaShear    !< If true, use the kappa_shear module to find the
                              !! shear-driven diapycnal diffusivity.
  logical :: Vertex_Shear     !< If true, do the calculations of the shear-driven mixing
                              !! at the cell vertices (i.e., the vorticity points).
  logical :: use_CVMix_shear  !< If true, use one of the CVMix modules to find
                              !! shear-driven diapycnal diffusivity.
  logical :: double_diffusion !< If true, enable double-diffusive mixing using an old method.
  logical :: use_CVMix_ddiff  !< If true, enable double-diffusive mixing via CVMix.
  logical :: use_tidal_mixing !< If true, activate tidal mixing diffusivity.
  logical :: simple_TKE_to_Kd !< If true, uses a simple estimate of Kd/TKE that
                              !! does not rely on a layer-formulation.
  real    :: Max_Rrho_salt_fingers      !< max density ratio for salt fingering
  real    :: Max_salt_diff_salt_fingers !< max salt diffusivity for salt fingers [Z2 T-1 ~> m2 s-1]
  real    :: Kv_molecular               !< molecular visc for double diff convect [Z2 T-1 ~> m2 s-1]

  logical :: answers_2018   !< If true, use the order of arithmetic and expressions that recover the
                            !! answers from the end of 2018.  Otherwise, use updated and more robust
                            !! forms of the same expressions.

  character(len=200) :: inputdir !< The directory in which input files are found
  type(user_change_diff_CS), pointer :: user_change_diff_CSp => NULL() !< Control structure for a child module
  type(Kappa_shear_CS),      pointer :: kappaShear_CSp       => NULL() !< Control structure for a child module
  type(CVMix_shear_cs),      pointer :: CVMix_shear_csp      => NULL() !< Control structure for a child module
  type(CVMix_ddiff_cs),      pointer :: CVMix_ddiff_csp      => NULL() !< Control structure for a child module
  type(bkgnd_mixing_cs),     pointer :: bkgnd_mixing_csp     => NULL() !< Control structure for a child module
  type(int_tide_CS),         pointer :: int_tide_CSp         => NULL() !< Control structure for a child module
  type(tidal_mixing_cs),     pointer :: tidal_mixing_CSp     => NULL() !< Control structure for a child module

  !>@{ Diagnostic IDs
  integer :: id_maxTKE     = -1, id_TKE_to_Kd   = -1, id_Kd_user    = -1
  integer :: id_Kd_layer   = -1, id_Kd_BBL      = -1, id_N2         = -1
  integer :: id_Kd_Work    = -1, id_KT_extra    = -1, id_KS_extra   = -1, id_R_rho = -1
  integer :: id_Kd_bkgnd   = -1, id_Kv_bkgnd    = -1
  !>@}

end type set_diffusivity_CS

!> This structure has memory for used in calculating diagnostics of diffusivity
type diffusivity_diags
  real, pointer, dimension(:,:,:) :: &
    N2_3d    => NULL(), & !< squared buoyancy frequency at interfaces [T-2 ~> s-2]
    Kd_user  => NULL(), & !< user-added diffusivity at interfaces [Z2 T-1 ~> m2 s-1]
    Kd_BBL   => NULL(), & !< BBL diffusivity at interfaces [Z2 T-1 ~> m2 s-1]
    Kd_work  => NULL(), & !< layer integrated work by diapycnal mixing [R Z3 T-3 ~> W m-2]
    maxTKE   => NULL(), & !< energy required to entrain to h_max [Z3 T-3 ~> m3 s-3]
    Kd_bkgnd => NULL(), & !< Background diffusivity at interfaces [Z2 T-1 ~> m2 s-1]
    Kv_bkgnd => NULL(), & !< Viscosity from ackground diffusivity at interfaces [Z2 T-1 ~> m2 s-1]
    KT_extra => NULL(), & !< double diffusion diffusivity for temp [Z2 T-1 ~> m2 s-1].
    KS_extra => NULL(), & !< double diffusion diffusivity for saln [Z2 T-1 ~> m2 s-1].
    drho_rat => NULL()    !< The density difference ratio used in double diffusion [nondim].
  real, pointer, dimension(:,:,:) :: TKE_to_Kd => NULL()
                          !< conversion rate (~1.0 / (G_Earth + dRho_lay)) between TKE
                          !! dissipated within a layer and Kd in that layer
                          !! [Z2 T-1 / Z3 T-3 = T2 Z-1 ~> s2 m-1]

end type diffusivity_diags

!>@{ CPU time clocks
integer :: id_clock_kappaShear, id_clock_CVMix_ddiff
!>@}

contains

!> Sets the interior vertical diffusion of scalars due to the following processes:
!! 1. Shear-driven mixing: two options, Jackson et at. and KPP interior;
!! 2. Background mixing via CVMix (Bryan-Lewis profile) or the scheme described by
!!    Harrison & Hallberg, JPO 2008;
!! 3. Double-diffusion, old method and new method via CVMix;
!! 4. Tidal mixing: many options available, see MOM_tidal_mixing.F90;
!! In addition, this subroutine has the option to set the interior vertical
!! viscosity associated with processes 1,2 and 4 listed above, which is stored in
!! visc%Kv_slow. Vertical viscosity due to shear-driven mixing is passed via
!! visc%Kv_shear
subroutine set_diffusivity(u, v, h, u_h, v_h, tv, fluxes, optics, visc, dt, &
                           G, GV, US, CS, Kd_lay, Kd_int, Kd_extra_T, Kd_extra_S)
  type(ocean_grid_type),     intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),     intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                             intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                             intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(in)    :: u_h  !< Zonal velocity interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(in)    :: v_h  !< Meridional velocity interpolated to h points [L T-1 ~> m s-1].
  type(thermo_var_ptrs),     intent(inout) :: tv   !< Structure with pointers to thermodynamic
                                                   !! fields. Out is for tv%TempxPmE.
  type(forcing),             intent(in)    :: fluxes !< A structure of thermodynamic surface fluxes
  type(optics_type),         pointer       :: optics !< A structure describing the optical
                                                   !!  properties of the ocean.
  type(vertvisc_type),       intent(inout) :: visc !< Structure containing vertical viscosities, bottom
                                                   !! boundary layer properies, and related fields.
  real,                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(set_diffusivity_CS),  pointer       :: CS   !< Module control structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                   optional, intent(out)   :: Kd_lay !< Diapycnal diffusivity of each layer [Z2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                   optional, intent(out)   :: Kd_int !< Diapycnal diffusivity at each interface [Z2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                   optional, intent(out)   :: Kd_extra_T !< The extra diffusivity at interfaces of
                                                     !! temperature due to double diffusion relative to
                                                     !! the diffusivity of density [Z2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                   optional, intent(out)   :: Kd_extra_S !< The extra diffusivity at interfaces of
                                                     !! salinity due to double diffusion relative to
                                                     !! the diffusivity of density [Z2 T-1 ~> m2 s-1].

  ! local variables
  real, dimension(SZI_(G)) :: &
    N2_bot        ! bottom squared buoyancy frequency [T-2 ~> s-2]

  type(diffusivity_diags)  :: dd ! structure with arrays of available diags

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    T_f, S_f      ! Temperature and salinity [degC] and [ppt] with properties in massless layers
                  ! filled vertically by diffusion or the properties after full convective adjustment.

  real, dimension(SZI_(G),SZK_(GV)) :: &
    N2_lay, &     !< Squared buoyancy frequency associated with layers [T-2 ~> s-2]
    Kd_lay_2d, &  !< The layer diffusivities [Z2 T-1 ~> m2 s-1]
    maxTKE, &     !< Energy required to entrain to h_max [Z3 T-3 ~> m3 s-3]
    TKE_to_Kd     !< Conversion rate (~1.0 / (G_Earth + dRho_lay)) between
                  !< TKE dissipated within a layer and Kd in that layer
                  !< [Z2 T-1 / Z3 T-3 = T2 Z-1 ~> s2 m-1]

  real, dimension(SZI_(G),SZK_(GV)+1) :: &
    N2_int,   &   !< squared buoyancy frequency associated at interfaces [T-2 ~> s-2]
    Kd_int_2d, &  !< The interface diffusivities [Z2 T-1 ~> m2 s-1]
    Kv_bkgnd, &   !< The background diffusion related interface viscosities [Z2 T-1 ~> m2 s-1]
    dRho_int, &   !< Locally referenced potential density difference across interfaces [R ~> kg m-3]
    KT_extra, &   !< Double difusion diffusivity of temperature [Z2 T-1 ~> m2 s-1]
    KS_extra      !< Double difusion diffusivity of salinity [Z2 T-1 ~> m2 s-1]

  real :: dissip        ! local variable for dissipation calculations [Z2 R T-3 ~> W m-3]
  real :: Omega2        ! squared absolute rotation rate [T-2 ~> s-2]

  logical   :: use_EOS      ! If true, compute density from T/S using equation of state.
  logical   :: TKE_to_Kd_used ! If true, TKE_to_Kd and maxTKE need to be calculated.
  integer   :: kb(SZI_(G))  ! The index of the lightest layer denser than the
                            ! buffer layer, or -1 without a bulk mixed layer.
  logical   :: showCallTree ! If true, show the call tree.

  integer :: i, j, k, is, ie, js, je, nz, isd, ied, jsd, jed

  real      :: kappa_dt_fill ! diffusivity times a timestep used to fill massless layers [Z2 ~> m2]

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("set_diffusivity(), MOM_set_diffusivity.F90")

  if (.not.associated(CS)) call MOM_error(FATAL,"set_diffusivity: "//&
         "Module must be initialized before it is used.")

  if (CS%answers_2018) then
    ! These hard-coded dimensional parameters are being replaced.
    kappa_dt_fill = US%m_to_Z**2 * 1.e-3 * 7200.
  else
    kappa_dt_fill = CS%Kd_smooth * dt
  endif
  Omega2 = CS%omega * CS%omega

  use_EOS = associated(tv%eqn_of_state)

  if ((CS%use_CVMix_ddiff .or. CS%double_diffusion) .and. &
      .not.(present(Kd_extra_T) .and. present(Kd_extra_S))) &
    call MOM_error(FATAL, "set_diffusivity: both Kd_extra_T and Kd_extra_S must be present "//&
                          "when USE_CVMIX_DDIFF or DOUBLE_DIFFUSION are true.")

  TKE_to_Kd_used = (CS%use_tidal_mixing .or. CS%ML_radiation .or. &
                   (CS%bottomdraglaw .and. .not.CS%use_LOTW_BBL_diffusivity))

  ! Set Kd_lay, Kd_int and Kv_slow to constant values, mostly to fill the halos.
  if (present(Kd_lay)) Kd_lay(:,:,:) = CS%Kd
  if (present(Kd_int)) Kd_int(:,:,:) = CS%Kd
  if (present(Kd_extra_T)) Kd_extra_T(:,:,:) = 0.0
  if (present(Kd_extra_S)) Kd_extra_S(:,:,:) = 0.0
  if (associated(visc%Kv_slow)) visc%Kv_slow(:,:,:) = CS%Kv

  ! Set up arrays for diagnostics.

  if (CS%id_N2 > 0) then
    allocate(dd%N2_3d(isd:ied,jsd:jed,nz+1)) ; dd%N2_3d(:,:,:) = 0.0
  endif
  if (CS%id_Kd_user > 0) then
    allocate(dd%Kd_user(isd:ied,jsd:jed,nz+1)) ; dd%Kd_user(:,:,:) = 0.0
  endif
  if (CS%id_Kd_work > 0) then
    allocate(dd%Kd_work(isd:ied,jsd:jed,nz)) ; dd%Kd_work(:,:,:) = 0.0
  endif
  if (CS%id_maxTKE > 0) then
    allocate(dd%maxTKE(isd:ied,jsd:jed,nz)) ; dd%maxTKE(:,:,:) = 0.0
  endif
  if (CS%id_TKE_to_Kd > 0) then
    allocate(dd%TKE_to_Kd(isd:ied,jsd:jed,nz)) ; dd%TKE_to_Kd(:,:,:) = 0.0
  endif
  if ((CS%double_diffusion) .and. (CS%id_KT_extra > 0)) then
    allocate(dd%KT_extra(isd:ied,jsd:jed,nz+1)) ; dd%KT_extra(:,:,:) = 0.0
  endif
  if ((CS%double_diffusion) .and. (CS%id_KS_extra > 0)) then
    allocate(dd%KS_extra(isd:ied,jsd:jed,nz+1)) ; dd%KS_extra(:,:,:) = 0.0
  endif
  if (CS%id_R_rho > 0) then
    allocate(dd%drho_rat(isd:ied,jsd:jed,nz+1)) ; dd%drho_rat(:,:,:) = 0.0
  endif
  if (CS%id_Kd_BBL > 0) then
    allocate(dd%Kd_BBL(isd:ied,jsd:jed,nz+1)) ; dd%Kd_BBL(:,:,:) = 0.0
  endif

  if (CS%id_Kd_bkgnd > 0) then
    allocate(dd%Kd_bkgnd(isd:ied,jsd:jed,nz+1)) ; dd%Kd_bkgnd(:,:,:) = 0.
  endif
  if (CS%id_Kv_bkgnd > 0) then
    allocate(dd%Kv_bkgnd(isd:ied,jsd:jed,nz+1)) ; dd%Kv_bkgnd(:,:,:) = 0.
  endif

  ! set up arrays for tidal mixing diagnostics
  call setup_tidal_diagnostics(G, GV, CS%tidal_mixing_CSp)

  if (CS%useKappaShear) then
    if (CS%debug) then
      call hchksum_pair("before calc_KS [uv]_h", u_h, v_h, G%HI, scale=US%L_T_to_m_s)
    endif
    call cpu_clock_begin(id_clock_kappaShear)
    if (CS%Vertex_shear) then
      call full_convection(G, GV, US, h, tv, T_f, S_f, fluxes%p_surf, &
                           (GV%Z_to_H**2)*kappa_dt_fill, halo=1)

      call calc_kappa_shear_vertex(u, v, h, T_f, S_f, tv, fluxes%p_surf, visc%Kd_shear, &
                                   visc%TKE_turb, visc%Kv_shear_Bu, dt, G, GV, US, CS%kappaShear_CSp)
      if (associated(visc%Kv_shear)) visc%Kv_shear(:,:,:) = 0.0 ! needed for other parameterizations
      if (CS%debug) then
        call hchksum(visc%Kd_shear, "after calc_KS_vert visc%Kd_shear", G%HI, scale=US%Z2_T_to_m2_s)
        call Bchksum(visc%Kv_shear_Bu, "after calc_KS_vert visc%Kv_shear_Bu", G%HI, scale=US%Z2_T_to_m2_s)
        call Bchksum(visc%TKE_turb, "after calc_KS_vert visc%TKE_turb", G%HI, scale=US%Z_to_m**2*US%s_to_T**2)
      endif
    else
      ! Changes: visc%Kd_shear ;  Sets: visc%Kv_shear and visc%TKE_turb
      call calculate_kappa_shear(u_h, v_h, h, tv, fluxes%p_surf, visc%Kd_shear, visc%TKE_turb, &
                                 visc%Kv_shear, dt, G, GV, US, CS%kappaShear_CSp)
      if (CS%debug) then
        call hchksum(visc%Kd_shear, "after calc_KS visc%Kd_shear", G%HI, scale=US%Z2_T_to_m2_s)
        call hchksum(visc%Kv_shear, "after calc_KS visc%Kv_shear", G%HI, scale=US%Z2_T_to_m2_s)
        call hchksum(visc%TKE_turb, "after calc_KS visc%TKE_turb", G%HI, scale=US%Z_to_m**2*US%s_to_T**2)
      endif
    endif
    call cpu_clock_end(id_clock_kappaShear)
    if (showCallTree) call callTree_waypoint("done with calculate_kappa_shear (set_diffusivity)")
  elseif (CS%use_CVMix_shear) then
    !NOTE{BGR}: this needs to be cleaned up.  It works in 1D case, but has not been tested outside.
    call calculate_CVMix_shear(u_h, v_h, h, tv, visc%Kd_shear, visc%Kv_shear, G, GV, US, CS%CVMix_shear_CSp)
    if (CS%debug) then
      call hchksum(visc%Kd_shear, "after CVMix_shear visc%Kd_shear", G%HI, scale=US%Z2_T_to_m2_s)
      call hchksum(visc%Kv_shear, "after CVMix_shear visc%Kv_shear", G%HI, scale=US%Z2_T_to_m2_s)
    endif
  elseif (associated(visc%Kv_shear)) then
    visc%Kv_shear(:,:,:) = 0.0 ! needed if calculate_kappa_shear is not enabled
  endif

  ! Smooth the properties through massless layers.
  if (use_EOS) then
    if (CS%debug) then
      call hchksum(tv%T, "before vert_fill_TS tv%T",G%HI)
      call hchksum(tv%S, "before vert_fill_TS tv%S",G%HI)
      call hchksum(h, "before vert_fill_TS h",G%HI, scale=GV%H_to_m)
    endif
    call vert_fill_TS(h, tv%T, tv%S, kappa_dt_fill, T_f, S_f, G, GV, larger_h_denom=.true.)
    if (CS%debug) then
      call hchksum(tv%T, "after vert_fill_TS tv%T",G%HI)
      call hchksum(tv%S, "after vert_fill_TS tv%S",G%HI)
      call hchksum(h, "after vert_fill_TS h",G%HI, scale=GV%H_to_m)
    endif
  endif

  !   Calculate the diffusivities, Kd_lay and Kd_int, for each layer and interface.  This would
  ! be an appropriate place to add a depth-dependent parameterization or another explicit
  ! parameterization of Kd.

  !$OMP parallel do default(shared) private(dRho_int,N2_lay,Kd_lay_2d,Kd_int_2d,Kv_bkgnd,N2_int,&
  !$OMP                                     N2_bot,KT_extra,KS_extra,TKE_to_Kd,maxTKE,dissip,kb)
  do j=js,je

    ! Set up variables related to the stratification.
    call find_N2(h, tv, T_f, S_f, fluxes, j, G, GV, US, CS, dRho_int, N2_lay, N2_int, N2_bot)

    if (associated(dd%N2_3d)) then
      do K=1,nz+1 ; do i=is,ie ; dd%N2_3d(i,j,K) = N2_int(i,K) ; enddo ; enddo
    endif

    ! Add background mixing
    call calculate_bkgnd_mixing(h, tv, N2_lay, Kd_lay_2d, Kd_int_2d, Kv_bkgnd, j, G, GV, US, CS%bkgnd_mixing_csp)
    ! Update Kv and 3-d diffusivity diagnostics.
    if (associated(visc%Kv_slow)) then ; do K=1,nz+1 ; do i=is,ie
      visc%Kv_slow(i,j,K) = visc%Kv_slow(i,j,K) + Kv_bkgnd(i,K)
    enddo ; enddo ; endif
    if (CS%id_Kv_bkgnd > 0) then ; do K=1,nz+1 ; do i=is,ie
      dd%Kv_bkgnd(i,j,K) = Kv_bkgnd(i,K)
    enddo ; enddo ; endif
    if (CS%id_Kd_bkgnd > 0) then ; do K=1,nz+1 ; do i=is,ie
      dd%Kd_bkgnd(i,j,K) = Kd_int_2d(i,K)
    enddo ; enddo ; endif

    ! Double-diffusion (old method)
    if (CS%double_diffusion) then
      call double_diffusion(tv, h, T_f, S_f, j, G, GV, US, CS, KT_extra, KS_extra)
      ! One of Kd_extra_T and Kd_extra_S is always 0. Kd_extra_S is positive for salt fingering.
      ! Kd_extra_T is positive for double diffusive convection.
      do K=2,nz ; do i=is,ie
        if (KS_extra(i,K) > KT_extra(i,K)) then ! salt fingering
          Kd_lay_2d(i,k-1) = Kd_lay_2d(i,k-1) + 0.5 * KT_extra(i,K)
          Kd_lay_2d(i,k)   = Kd_lay_2d(i,k)   + 0.5 * KT_extra(i,K)
          Kd_extra_S(i,j,K) = (KS_extra(i,K) - KT_extra(i,K))
          Kd_extra_T(i,j,K) = 0.0
        elseif (KT_extra(i,K) > 0.0) then ! double-diffusive convection
          Kd_lay_2d(i,k-1) = Kd_lay_2d(i,k-1) + 0.5 * KS_extra(i,K)
          Kd_lay_2d(i,k)   = Kd_lay_2d(i,k)   + 0.5 * KS_extra(i,K)
          Kd_extra_T(i,j,K) = (KT_extra(i,K) - KS_extra(i,K))
          Kd_extra_S(i,j,K) = 0.0
        else ! There is no double diffusion at this interface.
          Kd_extra_T(i,j,K) = 0.0
          Kd_extra_S(i,j,K) = 0.0
        endif
      enddo ; enddo
      if (associated(dd%KT_extra)) then ; do K=1,nz+1 ; do i=is,ie
        dd%KT_extra(i,j,K) = KT_extra(i,K)
      enddo ; enddo ; endif

      if (associated(dd%KS_extra)) then ; do K=1,nz+1 ; do i=is,ie
        dd%KS_extra(i,j,K) = KS_extra(i,K)
      enddo ; enddo ; endif
    endif

    ! Apply double diffusion via CVMix
    ! GMM, we need to pass HBL to compute_ddiff_coeffs, but it is not yet available.
    if (CS%use_CVMix_ddiff) then
      call cpu_clock_begin(id_clock_CVMix_ddiff)
      if (associated(dd%drho_rat)) then
        call compute_ddiff_coeffs(h, tv, G, GV, US, j, Kd_extra_T, Kd_extra_S, &
                                  CS%CVMix_ddiff_csp, dd%drho_rat)
      else
        call compute_ddiff_coeffs(h, tv, G, GV, US, j, Kd_extra_T, Kd_extra_S, CS%CVMix_ddiff_csp)
      endif
      call cpu_clock_end(id_clock_CVMix_ddiff)
    endif

    ! Calculate conversion ratios from TKE to layer diffusivities.
    if (TKE_to_Kd_used) then
      call find_TKE_to_Kd(h, tv, dRho_int, N2_lay, j, dt, G, GV, US, CS, TKE_to_Kd, maxTKE, kb)
      if (associated(dd%maxTKE)) then ; do k=1,nz ; do i=is,ie
        dd%maxTKE(i,j,k) = maxTKE(i,k)
      enddo ; enddo ; endif
      if (associated(dd%TKE_to_Kd)) then ; do k=1,nz ; do i=is,ie
        dd%TKE_to_Kd(i,j,k) = TKE_to_Kd(i,k)
      enddo ; enddo ; endif
    endif

    ! Add the input turbulent diffusivity.
    if (CS%useKappaShear .or. CS%use_CVMix_shear) then
      if (present(Kd_int)) then
        do K=2,nz ; do i=is,ie
          Kd_int_2d(i,K) = visc%Kd_shear(i,j,K) + 0.5 * (Kd_lay_2d(i,k-1) + Kd_lay_2d(i,k))
        enddo ; enddo
        do i=is,ie
          Kd_int_2d(i,1) = visc%Kd_shear(i,j,1) ! This isn't actually used. It could be 0.
          Kd_int_2d(i,nz+1) = 0.0
        enddo
      endif
      do k=1,nz ; do i=is,ie
        Kd_lay_2d(i,k) = Kd_lay_2d(i,k) + 0.5 * (visc%Kd_shear(i,j,K) + visc%Kd_shear(i,j,K+1))
      enddo ; enddo
    else
      if (present(Kd_int)) then
        do i=is,ie
          Kd_int_2d(i,1) = Kd_lay_2d(i,1) ; Kd_int_2d(i,nz+1) = 0.0
        enddo
        do K=2,nz ; do i=is,ie
          Kd_int_2d(i,K) = 0.5 * (Kd_lay_2d(i,k-1) + Kd_lay_2d(i,k))
        enddo ; enddo
      endif
    endif

    if (present(Kd_int)) then
      ! Add the ML_Rad diffusivity.
      if (CS%ML_radiation) &
        call add_MLrad_diffusivity(h, fluxes, j, G, GV, US, CS, TKE_to_Kd, Kd_lay_2d, Kd_int_2d)

      ! Add the Nikurashin and / or tidal bottom-driven mixing
      if (CS%use_tidal_mixing) &
        call calculate_tidal_mixing(h, N2_bot, j, TKE_to_Kd, maxTKE, G, GV, US, CS%tidal_mixing_CSp, &
                                    N2_lay, N2_int, Kd_lay_2d, Kd_int_2d, CS%Kd_max, visc%Kv_slow)

      ! This adds the diffusion sustained by the energy extracted from the flow by the bottom drag.
      if (CS%bottomdraglaw .and. (CS%BBL_effic>0.0)) then
        if (CS%use_LOTW_BBL_diffusivity) then
          call add_LOTW_BBL_diffusivity(h, u, v, tv, fluxes, visc, j, N2_int, G, GV, US, CS,  &
                                        dd%Kd_BBL, Kd_lay_2d, Kd_int_2d)
        else
          call add_drag_diffusivity(h, u, v,  tv, fluxes, visc, j, TKE_to_Kd, &
                                    maxTKE, kb, G, GV, US, CS, Kd_lay_2d, Kd_int_2d, dd%Kd_BBL)
        endif
      endif

      if (CS%limit_dissipation) then
        ! This calculates the dissipation ONLY from Kd calculated in this routine
        ! dissip has units of W/m3 (= kg/m3 * m2/s * 1/s2)
        !   1) a global constant,
        !   2) a dissipation proportional to N (aka Gargett) and
        !   3) dissipation corresponding to a (nearly) constant diffusivity.
        do K=2,nz ; do i=is,ie
          dissip = max( CS%dissip_min, &   ! Const. floor on dissip.
                        CS%dissip_N0 + CS%dissip_N1 * sqrt(N2_int(i,K)), & ! Floor aka Gargett
                        CS%dissip_N2 * N2_int(i,K)) ! Floor of Kd_min*rho0/F_Ri
          Kd_int_2d(i,K) = max(Kd_int_2d(i,K) , &  ! Apply floor to Kd
                              dissip * (CS%FluxRi_max / (GV%Rho0 * (N2_int(i,K) + Omega2))))
        enddo ; enddo
      endif

      ! Optionally add a uniform diffusivity at the interfaces.
      if (CS%Kd_add > 0.0) then ; do K=1,nz+1 ; do i=is,ie
        Kd_int_2d(i,K) = Kd_int_2d(i,K) + CS%Kd_add
      enddo ; enddo ; endif

      ! Copy the 2-d slices into the 3-d array that is exported.
      do K=1,nz+1 ; do i=is,ie
        Kd_int(i,j,K) = Kd_int_2d(i,K)
      enddo ; enddo

    else ! Kd_int is not present.

      ! Add the ML_Rad diffusivity.
      if (CS%ML_radiation) &
        call add_MLrad_diffusivity(h, fluxes, j, G, GV, US, CS, TKE_to_Kd, Kd_lay_2d)

      ! Add the Nikurashin and / or tidal bottom-driven mixing
      if (CS%use_tidal_mixing) &
        call calculate_tidal_mixing(h, N2_bot, j, TKE_to_Kd, maxTKE, G, GV, US, CS%tidal_mixing_CSp, &
                                    N2_lay, N2_int, Kd_lay_2d, Kd_max=CS%Kd_max, Kv=visc%Kv_slow)

      ! This adds the diffusion sustained by the energy extracted from the flow by the bottom drag.
      if (CS%bottomdraglaw .and. (CS%BBL_effic>0.0)) then
        if (CS%use_LOTW_BBL_diffusivity) then
          call add_LOTW_BBL_diffusivity(h, u, v, tv, fluxes, visc, j, N2_int, G, GV, US, CS,  &
                                        dd%Kd_BBL, Kd_lay_2d)
        else
          call add_drag_diffusivity(h, u, v,  tv, fluxes, visc, j, TKE_to_Kd, &
                                    maxTKE, kb, G, GV, US, CS, Kd_lay_2d, Kd_BBL=dd%Kd_BBL)
        endif
      endif

    endif

    if (CS%limit_dissipation) then
      ! This calculates the layer dissipation ONLY from Kd calculated in this routine
      ! dissip has units of W/m3 (= kg/m3 * m2/s * 1/s2)
      !   1) a global constant,
      !   2) a dissipation proportional to N (aka Gargett) and
      !   3) dissipation corresponding to a (nearly) constant diffusivity.
      do k=2,nz-1 ; do i=is,ie
        dissip = max( CS%dissip_min, &   ! Const. floor on dissip.
                      CS%dissip_N0 + CS%dissip_N1 * sqrt(N2_lay(i,k)), & ! Floor aka Gargett
                      CS%dissip_N2 * N2_lay(i,k)) ! Floor of Kd_min*rho0/F_Ri
        Kd_lay_2d(i,k) = max(Kd_lay_2d(i,k) , &  ! Apply floor to Kd
                            dissip * (CS%FluxRi_max / (GV%Rho0 * (N2_lay(i,k) + Omega2))))
      enddo ; enddo
    endif

    if (associated(dd%Kd_work)) then
      do k=1,nz ; do i=is,ie
        dd%Kd_Work(i,j,k) = GV%Rho0 * Kd_lay_2d(i,k) * N2_lay(i,k) * &
                            GV%H_to_Z*h(i,j,k)  ! Watt m-2 s = kg s-3
      enddo ; enddo
    endif

    ! Optionally add a uniform diffusivity to the layers.
    if ((CS%Kd_add > 0.0) .and. (present(Kd_lay))) then
      do k=1,nz ; do i=is,ie
        Kd_lay_2d(i,k) = Kd_lay_2d(i,k) + CS%Kd_add
      enddo ; enddo
    endif

    ! Copy the 2-d slices into the 3-d array that is exported; this was done above for Kd_int.
    if (present(Kd_lay)) then ; do k=1,nz ; do i=is,ie
      Kd_lay(i,j,k) = Kd_lay_2d(i,k)
    enddo ; enddo ; endif
  enddo ! j-loop

  if (CS%user_change_diff) then
    call user_change_diff(h, tv, G, GV, US, CS%user_change_diff_CSp, Kd_lay, Kd_int, &
                          T_f, S_f, dd%Kd_user)
  endif

  if (CS%debug) then
    if (present(Kd_lay)) call hchksum(Kd_lay, "Kd_lay", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)

    if (CS%useKappaShear) call hchksum(visc%Kd_shear, "Turbulent Kd", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)

    if (CS%use_CVMix_ddiff) then
      call hchksum(Kd_extra_T, "MOM_set_diffusivity: Kd_extra_T", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
      call hchksum(Kd_extra_S, "MOM_set_diffusivity: Kd_extra_S", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
    endif

    if (associated(visc%kv_bbl_u) .and. associated(visc%kv_bbl_v)) then
      call uvchksum("BBL Kv_bbl_[uv]", visc%kv_bbl_u, visc%kv_bbl_v, G%HI, &
                    haloshift=0, symmetric=.true., scale=US%Z2_T_to_m2_s, &
                    scalar_pair=.true.)
    endif

    if (associated(visc%bbl_thick_u) .and. associated(visc%bbl_thick_v)) then
      call uvchksum("BBL bbl_thick_[uv]", visc%bbl_thick_u, visc%bbl_thick_v, &
                    G%HI, haloshift=0, symmetric=.true., scale=US%Z_to_m, &
                    scalar_pair=.true.)
    endif

    if (associated(visc%Ray_u) .and. associated(visc%Ray_v)) then
      call uvchksum("Ray_[uv]", visc%Ray_u, visc%Ray_v, G%HI, 0, symmetric=.true., scale=US%Z_to_m*US%s_to_T)
    endif

  endif

  ! post diagnostics
  if (present(Kd_lay) .and. (CS%id_Kd_layer > 0)) call post_data(CS%id_Kd_layer, Kd_lay, CS%diag)

  ! background mixing
  if (CS%id_Kd_bkgnd > 0) call post_data(CS%id_Kd_bkgnd, dd%Kd_bkgnd, CS%diag)
  if (CS%id_Kv_bkgnd > 0) call post_data(CS%id_Kv_bkgnd, dd%Kv_bkgnd, CS%diag)

  ! tidal mixing
  call post_tidal_diagnostics(G, GV, h, CS%tidal_mixing_CSp)
  if (CS%id_N2 > 0)         call post_data(CS%id_N2,        dd%N2_3d,     CS%diag)
  if (CS%id_Kd_Work > 0)    call post_data(CS%id_Kd_Work,   dd%Kd_Work,   CS%diag)
  if (CS%id_maxTKE > 0)     call post_data(CS%id_maxTKE,    dd%maxTKE,    CS%diag)
  if (CS%id_TKE_to_Kd > 0)  call post_data(CS%id_TKE_to_Kd, dd%TKE_to_Kd, CS%diag)

  if (CS%id_Kd_user > 0)    call post_data(CS%id_Kd_user,   dd%Kd_user,   CS%diag)

  ! double diffusive mixing
  if (CS%double_diffusion) then
    if (CS%id_KT_extra > 0) call post_data(CS%id_KT_extra, dd%KT_extra, CS%diag)
    if (CS%id_KS_extra > 0) call post_data(CS%id_KS_extra, dd%KS_extra, CS%diag)
  elseif (CS%use_CVMix_ddiff) then
    if (CS%id_KT_extra > 0) call post_data(CS%id_KT_extra, Kd_extra_T, CS%diag)
    if (CS%id_KS_extra > 0) call post_data(CS%id_KS_extra, Kd_extra_S, CS%diag)
    if (CS%id_R_rho > 0) call post_data(CS%id_R_rho, dd%drho_rat, CS%diag)
  endif
  if (CS%id_Kd_BBL > 0)   call post_data(CS%id_Kd_BBL, dd%Kd_BBL, CS%diag)

  if (associated(dd%N2_3d)) deallocate(dd%N2_3d)
  if (associated(dd%Kd_work)) deallocate(dd%Kd_work)
  if (associated(dd%Kd_user)) deallocate(dd%Kd_user)
  if (associated(dd%maxTKE)) deallocate(dd%maxTKE)
  if (associated(dd%TKE_to_Kd)) deallocate(dd%TKE_to_Kd)
  if (associated(dd%KT_extra)) deallocate(dd%KT_extra)
  if (associated(dd%KS_extra)) deallocate(dd%KS_extra)
  if (associated(dd%drho_rat)) deallocate(dd%drho_rat)
  if (associated(dd%Kd_BBL)) deallocate(dd%Kd_BBL)

  if (showCallTree) call callTree_leave("set_diffusivity()")

end subroutine set_diffusivity

!> Convert turbulent kinetic energy to diffusivity
subroutine find_TKE_to_Kd(h, tv, dRho_int, N2_lay, j, dt, G, GV, US, CS, &
                          TKE_to_Kd, maxTKE, kb)
  type(ocean_grid_type),            intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),          intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),            intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),            intent(in)    :: tv   !< Structure containing pointers to any available
                                                          !! thermodynamic fields.
  real, dimension(SZI_(G),SZK_(GV)+1), intent(in) :: dRho_int !< Change in locally referenced potential density
                                                          !! across each interface [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: N2_lay !< The squared buoyancy frequency of the
                                                          !! layers [T-2 ~> s-2].
  integer,                          intent(in)    :: j    !< j-index of row to work on
  real,                             intent(in)    :: dt   !< Time increment [T ~> s].
  type(set_diffusivity_CS),         pointer       :: CS   !< Diffusivity control structure
  real, dimension(SZI_(G),SZK_(GV)), intent(out)  :: TKE_to_Kd !< The conversion rate between the
                                                          !! TKE dissipated within a layer and the
                                                          !! diapycnal diffusivity witin that layer,
                                                          !! usually (~Rho_0 / (G_Earth * dRho_lay))
                                                          !! [Z2 T-1 / Z3 T-3 = T2 Z-1 ~> s2 m-1]
  real, dimension(SZI_(G),SZK_(GV)), intent(out)  :: maxTKE !< The energy required to for a layer to entrain
                                                          !! to its maximum realizable thickness [Z3 T-3 ~> m3 s-3]
  integer, dimension(SZI_(G)),      intent(out)   :: kb   !< Index of lightest layer denser than the buffer
                                                          !! layer, or -1 without a bulk mixed layer.
  ! Local variables
  real, dimension(SZI_(G),SZK_(GV)) :: &
    ds_dsp1, &    ! coordinate variable (sigma-2) difference across an
                  ! interface divided by the difference across the interface
                  ! below it [nondim]
    dsp1_ds, &    ! inverse coordinate variable (sigma-2) difference
                  ! across an interface times the difference across the
                  ! interface above it [nondim]
    rho_0,   &    ! Layer potential densities relative to surface pressure [R ~> kg m-3]
    maxEnt        ! maxEnt is the maximum value of entrainment from below (with
                  ! compensating entrainment from above to keep the layer
                  ! density from changing) that will not deplete all of the
                  ! layers above or below a layer within a timestep [Z ~> m].
  real, dimension(SZI_(G)) :: &
    htot,    &    ! total thickness above or below a layer, or the
                  ! integrated thickness in the BBL [Z ~> m].
    mFkb,    &    ! total thickness in the mixed and buffer layers times ds_dsp1 [Z ~> m].
    p_ref,   &    ! array of tv%P_Ref pressures [R L2 T-2 ~> Pa]
    Rcv_kmb, &    ! coordinate density in the lowest buffer layer [R ~> kg m-3]
    p_0           ! An array of 0 pressures [R L2 T-2 ~> Pa]

  real :: dh_max      ! maximum amount of entrainment a layer could
                      ! undergo before entraining all fluid in the layers
                      ! above or below [Z ~> m].
  real :: dRho_lay    ! density change across a layer [R ~> kg m-3]
  real :: Omega2      ! rotation rate squared [T-2 ~> s-2]
  real :: G_Rho0      ! gravitation accel divided by Bouss ref density [Z T-2 R-1 ~> m4 s-2 kg-1]
  real :: G_IRho0     ! Alternate calculation of G_Rho0 for reproducibility [Z T-2 R-1 ~> m4 s-2 kg-1]
  real :: I_Rho0      ! inverse of Boussinesq reference density [R-1 ~> m3 kg-1]
  real :: I_dt        ! 1/dt [T-1 ~> s-1]
  real :: H_neglect   ! negligibly small thickness [H ~> m or kg m-2]
  real :: hN2pO2      ! h (N^2 + Omega^2), in [Z T-2 ~> m s-2].
  logical :: do_i(SZI_(G))

  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, k, is, ie, nz, i_rem, kmb, kb_min
  is = G%isc ; ie = G%iec ; nz = GV%ke

  I_dt      = 1.0 / dt
  Omega2    = CS%omega**2
  H_neglect = GV%H_subroundoff
  G_Rho0    = (US%L_to_Z**2 * GV%g_Earth) / (GV%Rho0)
  if (CS%answers_2018) then
    I_Rho0    = 1.0 / (GV%Rho0)
    G_IRho0 = (US%L_to_Z**2 * GV%g_Earth) * I_Rho0
  else
    G_IRho0 = G_Rho0
  endif

  ! Simple but coordinate-independent estimate of Kd/TKE
  if (CS%simple_TKE_to_Kd) then
    do k=1,nz ; do i=is,ie
      hN2pO2 = (GV%H_to_Z * h(i,j,k)) * (N2_lay(i,k) + Omega2) ! Units of Z T-2.
      if (hN2pO2>0.) then
        TKE_to_Kd(i,k) = 1.0 / hN2pO2 ! Units of T2 Z-1.
      else; TKE_to_Kd(i,k) = 0.; endif
      ! The maximum TKE conversion we allow is really a statement
      ! about the upper diffusivity we allow. Kd_max must be set.
      maxTKE(i,k) = hN2pO2 * CS%Kd_max ! Units of Z3 T-3.
    enddo ; enddo
    kb(is:ie) = -1 ! kb should not be used by any code in non-layered mode -AJA
    return
  endif

  ! Determine kb - the index of the shallowest active interior layer.
  if (CS%bulkmixedlayer) then
    kmb = GV%nk_rho_varies
    do i=is,ie ; p_0(i) = 0.0 ; p_ref(i) = tv%P_Ref ; enddo
    EOSdom(:) = EOS_domain(G%HI)
    do k=1,nz
      call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_0, rho_0(:,k), tv%eqn_of_state, EOSdom)
    enddo
    call calculate_density(tv%T(:,j,kmb), tv%S(:,j,kmb), p_ref, Rcv_kmb, tv%eqn_of_state, EOSdom)

    kb_min = kmb+1
    do i=is,ie
      !   Determine the next denser layer than the buffer layer in the
      ! coordinate density (sigma-2).
      do k=kmb+1,nz-1 ; if (Rcv_kmb(i) <= GV%Rlay(k)) exit ; enddo
      kb(i) = k

    !   Backtrack, in case there are massive layers above that are stable
    ! in sigma-0.
      do k=kb(i)-1,kmb+1,-1
        if (rho_0(i,kmb) > rho_0(i,k)) exit
        if (h(i,j,k)>2.0*GV%Angstrom_H) kb(i) = k
      enddo
    enddo

    call set_density_ratios(h, tv, kb, G, GV, US, CS, j, ds_dsp1, rho_0)
  else ! not bulkmixedlayer
    kb_min = 2 ; kmb = 0
    do i=is,ie ; kb(i) = 1 ; enddo
    call set_density_ratios(h, tv, kb, G, GV, US, CS, j, ds_dsp1)
  endif

  ! Determine maxEnt - the maximum permitted entrainment from below by each
  ! interior layer.
  do k=2,nz-1 ; do i=is,ie
    dsp1_ds(i,k) = 1.0 / ds_dsp1(i,k)
  enddo ; enddo
  do i=is,ie ; dsp1_ds(i,nz) = 0.0 ; enddo

  if (CS%bulkmixedlayer) then
    kmb = GV%nk_rho_varies
    do i=is,ie
      htot(i) = GV%H_to_Z*h(i,j,kmb)
      mFkb(i) = 0.0
      if (kb(i) < nz) &
        mFkb(i) = ds_dsp1(i,kb(i)) * (GV%H_to_Z*(h(i,j,kmb) - GV%Angstrom_H))
    enddo
    do k=1,kmb-1 ; do i=is,ie
      htot(i) = htot(i) + GV%H_to_Z*h(i,j,k)
      mFkb(i) = mFkb(i) + ds_dsp1(i,k+1)*(GV%H_to_Z*(h(i,j,k) - GV%Angstrom_H))
    enddo ; enddo
  else
    do i=is,i
      maxEnt(i,1) = 0.0 ; htot(i) = GV%H_to_Z*(h(i,j,1) - GV%Angstrom_H)
    enddo
  endif
  do k=kb_min,nz-1 ; do i=is,ie
    if (k == kb(i)) then
      maxEnt(i,kb(i)) = mFkb(i)
    elseif (k > kb(i)) then
      if (CS%answers_2018) then
        maxEnt(i,k) = (1.0/dsp1_ds(i,k))*(maxEnt(i,k-1) + htot(i))
      else
        maxEnt(i,k) = ds_dsp1(i,k)*(maxEnt(i,k-1) + htot(i))
      endif
      htot(i) = htot(i) + GV%H_to_Z*(h(i,j,k) - GV%Angstrom_H)
    endif
  enddo ; enddo

  do i=is,ie
    htot(i) = GV%H_to_Z*(h(i,j,nz) - GV%Angstrom_H) ; maxEnt(i,nz) = 0.0
    do_i(i) = (G%mask2dT(i,j) > 0.5)
  enddo
  do k=nz-1,kb_min,-1
    i_rem = 0
    do i=is,ie ; if (do_i(i)) then
      if (k<kb(i)) then ; do_i(i) = .false. ; cycle ; endif
      i_rem = i_rem + 1  ! Count the i-rows that are still being worked on.
      maxEnt(i,k) = MIN(maxEnt(i,k),dsp1_ds(i,k+1)*maxEnt(i,k+1) + htot(i))
      htot(i) = htot(i) + GV%H_to_Z*(h(i,j,k) - GV%Angstrom_H)
    endif ; enddo
    if (i_rem == 0) exit
  enddo ! k-loop

  ! Now set maxTKE and TKE_to_Kd.
  do i=is,ie
    maxTKE(i,1) = 0.0 ; TKE_to_Kd(i,1) = 0.0
    maxTKE(i,nz) = 0.0 ; TKE_to_Kd(i,nz) = 0.0
  enddo
  do k=2,kmb ; do i=is,ie
    maxTKE(i,k) = 0.0
    TKE_to_Kd(i,k) = 1.0 / ((N2_lay(i,k) + Omega2) * &
                            (GV%H_to_Z*(h(i,j,k) + H_neglect)))
  enddo ; enddo
  do k=kmb+1,kb_min-1 ; do i=is,ie
    !   These are the properties in the deeper mixed and buffer layers, and
    ! should perhaps be revisited.
    maxTKE(i,k) = 0.0 ; TKE_to_Kd(i,k) = 0.0
  enddo ; enddo
  do k=kb_min,nz-1 ; do i=is,ie
    if (k<kb(i)) then
      maxTKE(i,k) = 0.0
      TKE_to_Kd(i,k) = 0.0
    else
      ! maxTKE is found by determining the kappa that gives maxEnt.
      !  kappa_max = I_dt * dRho_int(i,K+1) * maxEnt(i,k) * &
      !             (GV%H_to_Z*h(i,j,k) + dh_max) / dRho_lay
      !  maxTKE(i,k) = (GV%g_Earth*US%L_to_Z**2) * dRho_lay * kappa_max
      ! dRho_int should already be non-negative, so the max is redundant?
      dh_max = maxEnt(i,k) * (1.0 + dsp1_ds(i,k))
      dRho_lay = 0.5 * max(dRho_int(i,K) + dRho_int(i,K+1), 0.0)
      maxTKE(i,k) = I_dt * (G_IRho0 * &
          (0.5*max(dRho_int(i,K+1) + dsp1_ds(i,k)*dRho_int(i,K), 0.0))) * &
           ((GV%H_to_Z*h(i,j,k) + dh_max) * maxEnt(i,k))
      ! TKE_to_Kd should be rho_InSitu / G_Earth * (delta rho_InSitu)
      ! The omega^2 term in TKE_to_Kd is due to a rescaling of the efficiency of turbulent
      ! mixing by a factor of N^2 / (N^2 + Omega^2), as proposed by Melet et al., 2013?
      TKE_to_Kd(i,k) = 1.0 / (G_Rho0 * dRho_lay + &
                              CS%omega**2 * GV%H_to_Z*(h(i,j,k) + H_neglect))
    endif
  enddo ; enddo

end subroutine find_TKE_to_Kd

!> Calculate Brunt-Vaisala frequency, N^2.
subroutine find_N2(h, tv, T_f, S_f, fluxes, j, G, GV, US, CS, dRho_int, &
                   N2_lay, N2_int, N2_bot)
  type(ocean_grid_type),    intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),  intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),    intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),    intent(in)  :: tv   !< Structure containing pointers to any available
                                                !! thermodynamic fields.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: T_f  !< layer temperature with the values in massless layers
                                                !! filled vertically by diffusion [degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: S_f  !< Layer salinities with values in massless
                                                !! layers filled vertically by diffusion [ppt].
  type(forcing),            intent(in)  :: fluxes !< A structure of thermodynamic surface fluxes
  integer,                  intent(in)  :: j    !< j-index of row to work on
  type(set_diffusivity_CS), pointer     :: CS   !< Diffusivity control structure
  real, dimension(SZI_(G),SZK_(GV)+1), &
                            intent(out) :: dRho_int !< Change in locally referenced potential density
                                                !! across each interface [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)+1), &
                            intent(out) :: N2_int !< The squared buoyancy frequency at the interfaces [T-2 ~> s-2].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(out) :: N2_lay !< The squared buoyancy frequency of the layers [T-2 ~> s-2].
  real, dimension(SZI_(G)), intent(out) :: N2_bot !< The near-bottom squared buoyancy frequency [T-2 ~> s-2].
  ! Local variables
  real, dimension(SZI_(G),SZK_(GV)+1) :: &
    dRho_int_unfilt, & ! unfiltered density differences across interfaces [R ~> kg m-3]
    dRho_dT,         & ! partial derivative of density wrt temp [R degC-1 ~> kg m-3 degC-1]
    dRho_dS            ! partial derivative of density wrt saln [R ppt-1 ~> kg m-3 ppt-1]

  real, dimension(SZI_(G)) :: &
    pres,      &  ! pressure at each interface [R L2 T-2 ~> Pa]
    Temp_int,  &  ! temperature at each interface [degC]
    Salin_int, &  ! salinity at each interface [ppt]
    drho_bot,  &  ! A density difference [R ~> kg m-3]
    h_amp,     &  ! The topographic roughness amplitude [Z ~> m].
    hb,        &  ! The thickness of the bottom layer [Z ~> m].
    z_from_bot    ! The hieght above the bottom [Z ~> m].

  real :: Rml_base  ! density of the deepest variable density layer
  real :: dz_int    ! thickness associated with an interface [Z ~> m].
  real :: G_Rho0    ! gravitation acceleration divided by Bouss reference density
                    ! times some unit conversion factors [Z T-2 R-1 ~> m4 s-2 kg-1].
  real :: H_neglect ! negligibly small thickness, in the same units as h.

  logical :: do_i(SZI_(G)), do_any
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, k, is, ie, nz

  is = G%isc ; ie = G%iec ; nz = GV%ke
  G_Rho0    = (US%L_to_Z**2 * GV%g_Earth) / (GV%Rho0)
  H_neglect = GV%H_subroundoff

  ! Find the (limited) density jump across each interface.
  do i=is,ie
    dRho_int(i,1) = 0.0 ; dRho_int(i,nz+1) = 0.0
    dRho_int_unfilt(i,1) = 0.0 ; dRho_int_unfilt(i,nz+1) = 0.0
  enddo
  if (associated(tv%eqn_of_state)) then
    if (associated(fluxes%p_surf)) then
      do i=is,ie ; pres(i) = fluxes%p_surf(i,j) ; enddo
    else
      do i=is,ie ; pres(i) = 0.0 ; enddo
    endif
    EOSdom(:) = EOS_domain(G%HI)
    do K=2,nz
      do i=is,ie
        pres(i) = pres(i) + (GV%g_Earth*GV%H_to_RZ)*h(i,j,k-1)
        Temp_Int(i) = 0.5 * (T_f(i,j,k) + T_f(i,j,k-1))
        Salin_Int(i) = 0.5 * (S_f(i,j,k) + S_f(i,j,k-1))
      enddo
      call calculate_density_derivs(Temp_int, Salin_int, pres, dRho_dT(:,K), dRho_dS(:,K), &
                                    tv%eqn_of_state, EOSdom)
      do i=is,ie
        dRho_int(i,K) = max(dRho_dT(i,K)*(T_f(i,j,k) - T_f(i,j,k-1)) + &
                            dRho_dS(i,K)*(S_f(i,j,k) - S_f(i,j,k-1)), 0.0)
        dRho_int_unfilt(i,K) = max(dRho_dT(i,K)*(tv%T(i,j,k) - tv%T(i,j,k-1)) + &
                            dRho_dS(i,K)*(tv%S(i,j,k) - tv%S(i,j,k-1)), 0.0)
      enddo
    enddo
  else
    do K=2,nz ; do i=is,ie
      dRho_int(i,K) = GV%Rlay(k) - GV%Rlay(k-1)
    enddo ; enddo
  endif

  ! Set the buoyancy frequencies.
  do k=1,nz ; do i=is,ie
    N2_lay(i,k) = G_Rho0 * 0.5*(dRho_int(i,K) + dRho_int(i,K+1)) / &
                  (GV%H_to_Z*(h(i,j,k) + H_neglect))
  enddo ; enddo
  do i=is,ie ; N2_int(i,1) = 0.0 ; N2_int(i,nz+1) = 0.0 ; enddo
  do K=2,nz ; do i=is,ie
    N2_int(i,K) = G_Rho0 * dRho_int(i,K) / &
                  (0.5*GV%H_to_Z*(h(i,j,k-1) + h(i,j,k) + H_neglect))
  enddo ; enddo

  ! Find the bottom boundary layer stratification, and use this in the deepest layers.
  do i=is,ie
    hb(i) = 0.0 ; dRho_bot(i) = 0.0 ; h_amp(i) = 0.0
    z_from_bot(i) = 0.5*GV%H_to_Z*h(i,j,nz)
    do_i(i) = (G%mask2dT(i,j) > 0.5)
  enddo
  if (CS%use_tidal_mixing) call tidal_mixing_h_amp(h_amp, G, j, CS%tidal_mixing_CSp)

  do k=nz,2,-1
    do_any = .false.
    do i=is,ie ; if (do_i(i)) then
      dz_int = 0.5*GV%H_to_Z*(h(i,j,k) + h(i,j,k-1))
      z_from_bot(i) = z_from_bot(i) + dz_int ! middle of the layer above

      hb(i) = hb(i) + dz_int
      drho_bot(i) = drho_bot(i) + dRho_int(i,K)

      if (z_from_bot(i) > h_amp(i)) then
        if (k>2) then
          ! Always include at least one full layer.
          hb(i) = hb(i) + 0.5*GV%H_to_Z*(h(i,j,k-1) + h(i,j,k-2))
          drho_bot(i) = drho_bot(i) + dRho_int(i,K-1)
        endif
        do_i(i) = .false.
      else
        do_any = .true.
      endif
    endif ; enddo
    if (.not.do_any) exit
  enddo

  do i=is,ie
    if (hb(i) > 0.0) then
      N2_bot(i) = (G_Rho0 * drho_bot(i)) / hb(i)
    else ;  N2_bot(i) = 0.0 ; endif
    z_from_bot(i) = 0.5*GV%H_to_Z*h(i,j,nz)
    do_i(i) = (G%mask2dT(i,j) > 0.5)
  enddo

  do k=nz,2,-1
    do_any = .false.
    do i=is,ie ; if (do_i(i)) then
      dz_int = 0.5*GV%H_to_Z*(h(i,j,k) + h(i,j,k-1))
      z_from_bot(i) = z_from_bot(i) + dz_int ! middle of the layer above

      N2_int(i,K) = N2_bot(i)
      if (k>2) N2_lay(i,k-1) = N2_bot(i)

      if (z_from_bot(i) > h_amp(i)) then
        if (k>2) N2_int(i,K-1) = N2_bot(i)
        do_i(i) = .false.
      else
        do_any = .true.
      endif
    endif ; enddo
    if (.not.do_any) exit
  enddo

  if (associated(tv%eqn_of_state)) then
    do K=1,nz+1 ; do i=is,ie
      dRho_int(i,K) = dRho_int_unfilt(i,K)
    enddo ; enddo
  endif

end subroutine find_N2

!> This subroutine sets the additional diffusivities of temperature and
!! salinity due to double diffusion, using the same functional form as is
!! used in MOM4.1, and taken from an NCAR technical note (REF?) that updates
!! what was in Large et al. (1994).  All the coefficients here should probably
!! be made run-time variables rather than hard-coded constants.
!!
!! \todo Find reference for NCAR tech note above.
subroutine double_diffusion(tv, h, T_f, S_f, j, G, GV, US, CS, Kd_T_dd, Kd_S_dd)
  type(ocean_grid_type),    intent(in)  :: G   !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)  :: US  !< A dimensional unit scaling type
  type(thermo_var_ptrs),    intent(in)  :: tv  !< Structure containing pointers to any available
                                               !! thermodynamic fields; absent fields have NULL ptrs.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: h   !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: T_f !< layer temperatures with the values in massless layers
                                               !! filled vertically by diffusion [degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: S_f !< Layer salinities with values in massless
                                               !! layers filled vertically by diffusion [ppt].
  integer,                  intent(in)  :: j   !< Meridional index upon which to work.
  type(set_diffusivity_CS), pointer     :: CS  !< Module control structure.
  real, dimension(SZI_(G),SZK_(GV)+1),       &
                            intent(out) :: Kd_T_dd !< Interface double diffusion diapycnal
                                               !! diffusivity for temp [Z2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1),       &
                            intent(out) :: Kd_S_dd !< Interface double diffusion diapycnal
                                               !! diffusivity for saln [Z2 T-1 ~> m2 s-1].

  real, dimension(SZI_(G)) :: &
    dRho_dT,  &    ! partial derivatives of density wrt temp [R degC-1 ~> kg m-3 degC-1]
    dRho_dS,  &    ! partial derivatives of density wrt saln [R ppt-1 ~> kg m-3 ppt-1]
    pres,     &    ! pressure at each interface [R L2 T-2 ~> Pa]
    Temp_int, &    ! temperature at interfaces [degC]
    Salin_int      ! Salinity at interfaces [ppt]

  real ::  alpha_dT ! density difference between layers due to temp diffs [R ~> kg m-3]
  real ::  beta_dS  ! density difference between layers due to saln diffs [R ~> kg m-3]

  real :: Rrho    ! vertical density ratio [nondim]
  real :: diff_dd ! factor for double-diffusion [nondim]
  real :: Kd_dd   ! The dominant double diffusive diffusivity [Z2 T-1 ~> m2 s-1]
  real :: prandtl ! flux ratio for diffusive convection regime

  real, parameter :: Rrho0  = 1.9 ! limit for double-diffusive density ratio [nondim]

  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, k, is, ie, nz
  is = G%isc ; ie = G%iec ; nz = GV%ke

  if (associated(tv%eqn_of_state)) then
    do i=is,ie
      pres(i) = 0.0 ; Kd_T_dd(i,1) = 0.0 ; Kd_S_dd(i,1) = 0.0
      Kd_T_dd(i,nz+1) = 0.0 ; Kd_S_dd(i,nz+1) = 0.0
    enddo
    if (associated(tv%p_surf)) then ; do i=is,ie ; pres(i) = tv%p_surf(i,j) ; enddo ; endif
    EOSdom(:) = EOS_domain(G%HI)
    do K=2,nz
      do i=is,ie
        pres(i) = pres(i) + (GV%g_Earth*GV%H_to_RZ)*h(i,j,k-1)
        Temp_Int(i) = 0.5 * (T_f(i,j,k-1) + T_f(i,j,k))
        Salin_Int(i) = 0.5 * (S_f(i,j,k-1) + S_f(i,j,k))
      enddo
      call calculate_density_derivs(Temp_int, Salin_int, pres, dRho_dT, dRho_dS, &
                                    tv%eqn_of_state, EOSdom)

      do i=is,ie
        alpha_dT = -1.0*dRho_dT(i) * (T_f(i,j,k-1) - T_f(i,j,k))
        beta_dS  = dRho_dS(i) * (S_f(i,j,k-1) - S_f(i,j,k))

        if ((alpha_dT > beta_dS) .and. (beta_dS > 0.0)) then  ! salt finger case
          Rrho = min(alpha_dT / beta_dS, Rrho0)
          diff_dd = 1.0 - ((RRho-1.0)/(RRho0-1.0))
          Kd_dd = CS%Max_salt_diff_salt_fingers * diff_dd*diff_dd*diff_dd
          Kd_T_dd(i,K) = 0.7 * Kd_dd
          Kd_S_dd(i,K) = Kd_dd
        elseif ((alpha_dT < 0.) .and. (beta_dS < 0.) .and. (alpha_dT > beta_dS)) then ! diffusive convection
          Rrho = alpha_dT / beta_dS
          Kd_dd = CS%Kv_molecular * 0.909 * exp(4.6 * exp(-0.54 * (1/Rrho - 1)))
          prandtl = 0.15*Rrho
          if (Rrho > 0.5) prandtl = (1.85-0.85/Rrho)*Rrho
          Kd_T_dd(i,K) = Kd_dd
          Kd_S_dd(i,K) = prandtl * Kd_dd
        else
          Kd_T_dd(i,K) = 0.0 ; Kd_S_dd(i,K) = 0.0
        endif
      enddo
    enddo
  endif

end subroutine double_diffusion

!> This routine adds diffusion sustained by flow energy extracted by bottom drag.
subroutine add_drag_diffusivity(h, u, v, tv, fluxes, visc, j, TKE_to_Kd, &
                                maxTKE, kb, G, GV, US, CS, Kd_lay, Kd_int, Kd_BBL)
  type(ocean_grid_type),            intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),          intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),            intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),            intent(in)    :: tv   !< Structure containing pointers to any available
                                                          !! thermodynamic fields.
  type(forcing),                    intent(in)    :: fluxes !< A structure of thermodynamic surface fluxes
  type(vertvisc_type),              intent(in)    :: visc !< Structure containing vertical viscosities, bottom
                                                          !! boundary layer properies, and related fields
  integer,                          intent(in)    :: j    !< j-index of row to work on
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: TKE_to_Kd !< The conversion rate between the TKE
                                                          !! TKE dissipated within  a layer and the
                                                          !! diapycnal diffusivity witin that layer,
                                                          !! usually (~Rho_0 / (G_Earth * dRho_lay))
                                                          !! [Z2 T-1 / Z3 T-3 = T2 Z-1 ~> s2 m-1]
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: maxTKE !< The energy required to for a layer to entrain
                                                          !! to its maximum-realizable thickness [Z3 T-3 ~> m3 s-3]
  integer, dimension(SZI_(G)),      intent(in)    :: kb   !< Index of lightest layer denser than the buffer
                                                          !! layer, or -1 without a bulk mixed layer
  type(set_diffusivity_CS),         pointer       :: CS   !< Diffusivity control structure
  real, dimension(SZI_(G),SZK_(GV)), intent(inout) :: Kd_lay !< The diapycnal diffusivity in layers,
                                                            !! [Z2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1), &
                          optional, intent(inout) :: Kd_int !< The diapycnal diffusivity at interfaces,
                                                            !! [Z2 T-1 ~> m2 s-1].
  real, dimension(:,:,:),           pointer       :: Kd_BBL !< Interface BBL diffusivity [Z2 T-1 ~> m2 s-1].

! This routine adds diffusion sustained by flow energy extracted by bottom drag.

  real, dimension(SZK_(GV)+1) :: &
    Rint          ! coordinate density of an interface [R ~> kg m-3]
  real, dimension(SZI_(G)) :: &
    htot, &       ! total thickness above or below a layer, or the
                  ! integrated thickness in the BBL [Z ~> m].
    rho_htot, &   ! running integral with depth of density [Z R ~> kg m-2]
    gh_sum_top, & ! BBL value of g'h that can be supported by
                  ! the local ustar, times R0_g [R ~> kg m-2]
    Rho_top, &    ! density at top of the BBL [R ~> kg m-3]
    TKE, &        ! turbulent kinetic energy available to drive
                  ! bottom-boundary layer mixing in a layer [Z3 T-3 ~> m3 s-3]
    I2decay       ! inverse of twice the TKE decay scale [Z-1 ~> m-1].

  real    :: TKE_to_layer   ! TKE used to drive mixing in a layer [Z3 T-3 ~> m3 s-3]
  real    :: TKE_Ray        ! TKE from layer Rayleigh drag used to drive mixing in layer [Z3 T-3 ~> m3 s-3]
  real    :: TKE_here       ! TKE that goes into mixing in this layer [Z3 T-3 ~> m3 s-3]
  real    :: dRl, dRbot     ! temporaries holding density differences [R ~> kg m-3]
  real    :: cdrag_sqrt     ! square root of the drag coefficient [nondim]
  real    :: ustar_h        ! value of ustar at a thickness point [Z T-1 ~> m s-1].
  real    :: absf           ! average absolute Coriolis parameter around a thickness point [T-1 ~> s-1]
  real    :: R0_g           ! Rho0 / G_Earth [R T2 Z-1 m-1 ~> kg s2 m-5]
  real    :: I_rho0         ! 1 / RHO0 [R-1 ~> m3 kg-1]
  real    :: delta_Kd       ! increment to Kd from the bottom boundary layer mixing [Z2 T-1 ~> m2 s-1].
  logical :: Rayleigh_drag  ! Set to true if Rayleigh drag velocities
                            ! defined in visc, on the assumption that this
                            ! extracted energy also drives diapycnal mixing.

  logical :: domore, do_i(SZI_(G))
  logical :: do_diag_Kd_BBL

  integer :: i, k, is, ie, nz, i_rem, kb_min
  is = G%isc ; ie = G%iec ; nz = GV%ke

  do_diag_Kd_BBL = associated(Kd_BBL)

  if (.not.(CS%bottomdraglaw .and. (CS%BBL_effic>0.0))) return

  cdrag_sqrt = sqrt(CS%cdrag)
  TKE_Ray = 0.0 ; Rayleigh_drag = .false.
  if (associated(visc%Ray_u) .and. associated(visc%Ray_v)) Rayleigh_drag = .true.

  I_Rho0 = 1.0 / (GV%Rho0)
  R0_g = GV%Rho0 / (US%L_to_Z**2 * GV%g_Earth)

  do K=2,nz ; Rint(K) = 0.5*(GV%Rlay(k-1)+GV%Rlay(k)) ; enddo

  kb_min = max(GV%nk_rho_varies+1,2)

  ! The turbulence decay scale is 0.5*ustar/f from K&E & MOM_vertvisc.F90
  ! Any turbulence that makes it into the mixed layers is assumed
  ! to be relatively small and is discarded.
  do i=is,ie
    ustar_h = visc%ustar_BBL(i,j)
    if (associated(fluxes%ustar_tidal)) &
      ustar_h = ustar_h + fluxes%ustar_tidal(i,j)
    absf = 0.25 * ((abs(G%CoriolisBu(I-1,J-1)) + abs(G%CoriolisBu(I,J))) + &
                   (abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J-1))))
    if ((ustar_h > 0.0) .and. (absf > 0.5*CS%IMax_decay*ustar_h))  then
      I2decay(i) = absf / ustar_h
    else
      ! The maximum decay scale should be something of order 200 m.
      ! If ustar_h = 0, this is land so this value doesn't matter.
      I2decay(i) = 0.5*CS%IMax_decay
    endif
    TKE(i) = ((CS%BBL_effic * cdrag_sqrt) * exp(-I2decay(i)*(GV%H_to_Z*h(i,j,nz))) ) * &
             visc%TKE_BBL(i,j)

    if (associated(fluxes%TKE_tidal)) &
      TKE(i) = TKE(i) + fluxes%TKE_tidal(i,j) * I_Rho0 * &
           (CS%BBL_effic * exp(-I2decay(i)*(GV%H_to_Z*h(i,j,nz))))

    ! Distribute the work over a BBL of depth 20^2 ustar^2 / g' following
    ! Killworth & Edwards (1999) and Zilitikevich & Mironov (1996).
    ! Rho_top is determined by finding the density where
    ! integral(bottom, Z) (rho(z') - rho(Z)) dz' = rho_0 400 ustar^2 / g

    gh_sum_top(i) = R0_g * 400.0 * ustar_h**2

    do_i(i) = (G%mask2dT(i,j) > 0.5)
    htot(i) = GV%H_to_Z*h(i,j,nz)
    rho_htot(i) = GV%Rlay(nz)*(GV%H_to_Z*h(i,j,nz))
    Rho_top(i) = GV%Rlay(1)
    if (CS%bulkmixedlayer .and. do_i(i)) Rho_top(i) = GV%Rlay(kb(i)-1)
  enddo

  do k=nz-1,2,-1 ; domore = .false.
    do i=is,ie ; if (do_i(i)) then
      htot(i) = htot(i) + GV%H_to_Z*h(i,j,k)
      rho_htot(i) = rho_htot(i) + GV%Rlay(k)*(GV%H_to_Z*h(i,j,k))
      if (htot(i)*GV%Rlay(k-1) <= (rho_htot(i) - gh_sum_top(i))) then
        ! The top of the mixing is in the interface atop the current layer.
        Rho_top(i) = (rho_htot(i) - gh_sum_top(i)) / htot(i)
        do_i(i) = .false.
      elseif (k <= kb(i)) then ; do_i(i) = .false.
      else ; domore = .true. ; endif
    endif ; enddo
    if (.not.domore) exit
  enddo ! k-loop

  do i=is,ie ; do_i(i) = (G%mask2dT(i,j) > 0.5) ; enddo
  do k=nz-1,kb_min,-1
    i_rem = 0
    do i=is,ie ; if (do_i(i)) then
      if (k<kb(i)) then ; do_i(i) = .false. ; cycle ; endif
      i_rem = i_rem + 1  ! Count the i-rows that are still being worked on.
      !   Apply vertical decay of the turbulent energy.  This energy is
      ! simply lost.
      TKE(i) = TKE(i) * exp(-I2decay(i) * (GV%H_to_Z*(h(i,j,k) + h(i,j,k+1))))

!      if (maxEnt(i,k) <= 0.0) cycle
      if (maxTKE(i,k) <= 0.0) cycle

  ! This is an analytic integral where diffusity is a quadratic function of
  ! rho that goes asymptotically to 0 at Rho_top (vaguely following KPP?).
      if (TKE(i) > 0.0) then
        if (Rint(K) <= Rho_top(i)) then
          TKE_to_layer = TKE(i)
        else
          dRl = Rint(K+1) - Rint(K) ; dRbot = Rint(K+1) - Rho_top(i)
          TKE_to_layer = TKE(i) * dRl * &
              (3.0*dRbot*(Rint(K) - Rho_top(i)) + dRl**2) / (dRbot**3)
        endif
      else ; TKE_to_layer = 0.0 ; endif

      ! TKE_Ray has been initialized to 0 above.
      if (Rayleigh_drag) TKE_Ray = 0.5*CS%BBL_effic * US%L_to_Z**2 * G%IareaT(i,j) * &
            ((G%areaCu(I-1,j) * visc%Ray_u(I-1,j,k) * u(I-1,j,k)**2 + &
              G%areaCu(I,j)   * visc%Ray_u(I,j,k)   * u(I,j,k)**2) + &
             (G%areaCv(i,J-1) * visc%Ray_v(i,J-1,k) * v(i,J-1,k)**2 + &
              G%areaCv(i,J)   * visc%Ray_v(i,J,k)   * v(i,J,k)**2))

      if (TKE_to_layer + TKE_Ray > 0.0) then
        if (CS%BBL_mixing_as_max) then
          if (TKE_to_layer + TKE_Ray > maxTKE(i,k)) &
              TKE_to_layer = maxTKE(i,k) - TKE_Ray

          TKE(i) = TKE(i) - TKE_to_layer

          if (Kd_lay(i,k) < (TKE_to_layer + TKE_Ray) * TKE_to_Kd(i,k)) then
            delta_Kd = (TKE_to_layer + TKE_Ray) * TKE_to_Kd(i,k) - Kd_lay(i,k)
            if ((CS%Kd_max >= 0.0) .and. (delta_Kd > CS%Kd_max)) then
              delta_Kd = CS%Kd_max
              Kd_lay(i,k) = Kd_lay(i,k) + delta_Kd
            else
              Kd_lay(i,k) = (TKE_to_layer + TKE_Ray) * TKE_to_Kd(i,k)
            endif
            if (present(Kd_int)) then
              Kd_int(i,K)   = Kd_int(i,K)   + 0.5 * delta_Kd
              Kd_int(i,K+1) = Kd_int(i,K+1) + 0.5 * delta_Kd
            endif
            if (do_diag_Kd_BBL) then
              Kd_BBL(i,j,K) = Kd_BBL(i,j,K) + 0.5 * delta_Kd
              Kd_BBL(i,j,K+1) = Kd_BBL(i,j,K+1) + 0.5 * delta_Kd
            endif
          endif
        else
          if (Kd_lay(i,k) >= maxTKE(i,k) * TKE_to_Kd(i,k)) then
            TKE_here = 0.0
            TKE(i) = TKE(i) + TKE_Ray
          elseif (Kd_lay(i,k) + (TKE_to_layer + TKE_Ray) * TKE_to_Kd(i,k) > &
                  maxTKE(i,k) * TKE_to_Kd(i,k)) then
            TKE_here = ((TKE_to_layer + TKE_Ray) + Kd_lay(i,k) / TKE_to_Kd(i,k)) - maxTKE(i,k)
            TKE(i) = (TKE(i) - TKE_here) + TKE_Ray
          else
            TKE_here = TKE_to_layer + TKE_Ray
            TKE(i) = TKE(i) - TKE_to_layer
          endif
          if (TKE(i) < 0.0) TKE(i) = 0.0 ! This should be unnecessary?

          if (TKE_here > 0.0) then
            delta_Kd = TKE_here * TKE_to_Kd(i,k)
            if (CS%Kd_max >= 0.0) delta_Kd = min(delta_Kd, CS%Kd_max)
            Kd_lay(i,k) = Kd_lay(i,k) + delta_Kd
            if (present(Kd_int)) then
              Kd_int(i,K)   = Kd_int(i,K)   + 0.5 * delta_Kd
              Kd_int(i,K+1) = Kd_int(i,K+1) + 0.5 * delta_Kd
            endif
            if (do_diag_Kd_BBL) then
              Kd_BBL(i,j,K) = Kd_BBL(i,j,K) + 0.5 * delta_Kd
              Kd_BBL(i,j,K+1) = Kd_BBL(i,j,K+1) + 0.5 * delta_Kd
            endif
          endif
        endif
      endif

      ! This may be risky - in the case that there are exactly zero
      ! velocities at 4 neighboring points, but nonzero velocities
      ! above the iterations would stop too soon. I don't see how this
      ! could happen in practice. RWH
      if ((TKE(i)<= 0.0) .and. (TKE_Ray == 0.0)) then
        do_i(i) = .false. ; i_rem = i_rem - 1
      endif

    endif ; enddo
    if (i_rem == 0) exit
  enddo ! k-loop

end subroutine add_drag_diffusivity

!> Calculates a BBL diffusivity use a Prandtl number 1 diffusivity with a law of the
!! wall turbulent viscosity, up to a BBL height where the energy used for mixing has
!! consumed the mechanical TKE input.
subroutine add_LOTW_BBL_diffusivity(h, u, v, tv, fluxes, visc, j, N2_int, &
                                    G, GV, US, CS, Kd_BBL, Kd_lay, Kd_int)
  type(ocean_grid_type),    intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),  intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),    intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: u  !< u component of flow [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(in)    :: v  !< v component of flow [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),    intent(in)    :: tv !< Structure containing pointers to any available
                                                !! thermodynamic fields.
  type(forcing),            intent(in)    :: fluxes !< Surface fluxes structure
  type(vertvisc_type),      intent(in)    :: visc !< Structure containing vertical viscosities, bottom
                                                  !! boundary layer properies, and related fields.
  integer,                  intent(in)    :: j  !< j-index of row to work on
  real, dimension(SZI_(G),SZK_(GV)+1), &
                            intent(in)    :: N2_int !< Square of Brunt-Vaisala at interfaces [T-2 ~> s-2]
  type(set_diffusivity_CS), pointer       :: CS !< Diffusivity control structure
  real, dimension(:,:,:),   pointer       :: Kd_BBL !< Interface BBL diffusivity [Z2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZK_(GV)), &
                  optional, intent(inout) :: Kd_lay !< Layer net diffusivity [Z2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZK_(GV)+1), &
                  optional, intent(inout) :: Kd_int !< Interface net diffusivity [Z2 T-1 ~> m2 s-1]

  ! Local variables
  real :: TKE_column       ! net TKE input into the column [Z3 T-3 ~> m3 s-3]
  real :: TKE_remaining    ! remaining TKE available for mixing in this layer and above [Z3 T-3 ~> m3 s-3]
  real :: TKE_consumed     ! TKE used for mixing in this layer [Z3 T-3 ~> m3 s-3]
  real :: TKE_Kd_wall      ! TKE associated with unlimited law of the wall mixing [Z3 T-3 ~> m3 s-3]
  real :: cdrag_sqrt       ! square root of the drag coefficient [nondim]
  real :: ustar            ! value of ustar at a thickness point [Z T-1 ~> m s-1].
  real :: ustar2           ! square of ustar, for convenience [Z2 T-2 ~> m2 s-2]
  real :: absf             ! average absolute value of Coriolis parameter around a thickness point [T-1 ~> s-1]
  real :: dh, dhm1         ! thickness of layers k and k-1, respecitvely [Z ~> m].
  real :: z_bot            ! distance to interface k from bottom [Z ~> m].
  real :: D_minus_z        ! distance to interface k from surface [Z ~> m].
  real :: total_thickness  ! total thickness of water column [Z ~> m].
  real :: Idecay           ! inverse of decay scale used for "Joule heating" loss of TKE with height [Z-1 ~> m-1].
  real :: Kd_wall          ! Law of the wall diffusivity [Z2 T-1 ~> m2 s-1].
  real :: Kd_lower         ! diffusivity for lower interface [Z2 T-1 ~> m2 s-1]
  real :: ustar_D          ! u* x D  [Z2 T-1 ~> m2 s-1].
  real :: I_Rho0           ! 1 / rho0 [R-1  ~> m3 kg-1]
  real :: N2_min           ! Minimum value of N2 to use in calculation of TKE_Kd_wall [T-2 ~> s-2]
  logical :: Rayleigh_drag ! Set to true if there are Rayleigh drag velocities defined in visc, on
                           ! the assumption that this extracted energy also drives diapycnal mixing.
  integer :: i, k, km1
  real, parameter :: von_karm = 0.41 ! Von Karman constant (http://en.wikipedia.org/wiki/Von_Karman_constant)
  logical :: do_diag_Kd_BBL

  if (.not.(CS%bottomdraglaw .and. (CS%BBL_effic>0.0))) return
  do_diag_Kd_BBL = associated(Kd_BBL)

  N2_min = 0.
  if (CS%LOTW_BBL_use_omega) N2_min = CS%omega**2

  ! Determine whether to add Rayleigh drag contribution to TKE
  Rayleigh_drag = .false.
  if (associated(visc%Ray_u) .and. associated(visc%Ray_v)) Rayleigh_drag = .true.
  I_Rho0 = 1.0 / (GV%Rho0)
  cdrag_sqrt = sqrt(CS%cdrag)

  do i=G%isc,G%iec ! Developed in single-column mode

    ! Column-wise parameters.
    absf = 0.25 * ((abs(G%CoriolisBu(I-1,J-1)) + abs(G%CoriolisBu(I,J))) + &
                   (abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J-1)))) ! Non-zero on equator!

    ! u* at the bottom [Z T-1 ~> m s-1].
    ustar = visc%ustar_BBL(i,j)
    ustar2 = ustar**2
    ! In add_drag_diffusivity(), fluxes%ustar_tidal is added in. This might be double counting
    ! since ustar_BBL should already include all contributions to u*? -AJA
    !### Examine the question of whether there is double counting of fluxes%ustar_tidal.
    if (associated(fluxes%ustar_tidal)) ustar = ustar + fluxes%ustar_tidal(i,j)

    ! The maximum decay scale should be something of order 200 m. We use the smaller of u*/f and
    ! (IMax_decay)^-1 as the decay scale. If ustar = 0, this is land so this value doesn't matter.
    Idecay = CS%IMax_decay
    if ((ustar > 0.0) .and. (absf > CS%IMax_decay * ustar)) Idecay = absf / ustar

    ! Energy input at the bottom [Z3 T-3 ~> m3 s-3].
    ! (Note that visc%TKE_BBL is in [Z3 T-3 ~> m3 s-3], set in set_BBL_TKE().)
    ! I am still unsure about sqrt(cdrag) in this expressions - AJA
    TKE_column = cdrag_sqrt * visc%TKE_BBL(i,j)
    ! Add in tidal dissipation energy at the bottom [R Z3 T-3 ~> m3 s-3].
    ! Note that TKE_tidal is in [R Z3 T-3 ~> W m-2].
    if (associated(fluxes%TKE_tidal)) &
      TKE_column = TKE_column + fluxes%TKE_tidal(i,j) * I_Rho0
    TKE_column = CS%BBL_effic * TKE_column ! Only use a fraction of the mechanical dissipation for mixing.

    TKE_remaining = TKE_column
    total_thickness = ( sum(h(i,j,:)) + GV%H_subroundoff )* GV%H_to_Z ! Total column thickness [Z ~> m].
    ustar_D = ustar * total_thickness
    z_bot = 0.
    Kd_lower = 0. ! Diffusivity on bottom boundary.

    ! Work upwards from the bottom, accumulating work used until it exceeds the available TKE input
    ! at the bottom.
    do k=GV%ke,2,-1
      dh = GV%H_to_Z * h(i,j,k) ! Thickness of this level [Z ~> m].
      km1 = max(k-1, 1)
      dhm1 = GV%H_to_Z * h(i,j,km1) ! Thickness of level above [Z ~> m].

      ! Add in additional energy input from bottom-drag against slopes (sides)
      if (Rayleigh_drag) TKE_remaining = TKE_remaining + &
            0.5*CS%BBL_effic * US%L_to_Z**2 * G%IareaT(i,j) * &
            ((G%areaCu(I-1,j) * visc%Ray_u(I-1,j,k) * u(I-1,j,k)**2 + &
              G%areaCu(I,j)   * visc%Ray_u(I,j,k)   * u(I,j,k)**2) + &
             (G%areaCv(i,J-1) * visc%Ray_v(i,J-1,k) * v(i,J-1,k)**2 + &
              G%areaCv(i,J)   * visc%Ray_v(i,J,k)   * v(i,J,k)**2))

      ! Exponentially decay TKE across the thickness of the layer.
      ! This is energy loss in addition to work done as mixing, apparently to Joule heating.
      TKE_remaining = exp(-Idecay*dh) * TKE_remaining

      z_bot = z_bot + h(i,j,k)*GV%H_to_Z ! Distance between upper interface of layer and the bottom [Z ~> m].
      D_minus_z = max(total_thickness - z_bot, 0.) ! Thickness above layer [Z ~> m].

      ! Diffusivity using law of the wall, limited by rotation, at height z [Z2 T-1 ~> m2 s-1].
      ! This calculation is at the upper interface of the layer
      if ( ustar_D + absf * ( z_bot * D_minus_z ) == 0.) then
        Kd_wall = 0.
      else
        Kd_wall = ((von_karm * ustar2) * (z_bot * D_minus_z)) &
                  / (ustar_D + absf * (z_bot * D_minus_z))
      endif

      ! TKE associated with Kd_wall [Z3 T-3 ~> m3 s-3].
      ! This calculation if for the volume spanning the interface.
      TKE_Kd_wall = Kd_wall * 0.5 * (dh + dhm1) * max(N2_int(i,k), N2_min)

      ! Now bound Kd such that the associated TKE is no greater than available TKE for mixing.
      if (TKE_Kd_wall > 0.) then
        TKE_consumed = min(TKE_Kd_wall, TKE_remaining)
        Kd_wall = (TKE_consumed / TKE_Kd_wall) * Kd_wall ! Scale Kd so that only TKE_consumed is used.
      else
        ! Either N2=0 or dh = 0.
        if (TKE_remaining > 0.) then
          Kd_wall = CS%Kd_max
        else
          Kd_wall = 0.
        endif
        TKE_consumed = 0.
      endif

      ! Now use up the appropriate about of TKE associated with the diffusivity chosen
      TKE_remaining = TKE_remaining - TKE_consumed ! Note this will be non-negative

      ! Add this BBL diffusivity to the model net diffusivity.
      if (present(Kd_int)) Kd_int(i,K) = Kd_int(i,K) + Kd_wall
      if (present(Kd_lay)) Kd_lay(i,k) = Kd_lay(i,k) + 0.5 * (Kd_wall + Kd_lower)
      Kd_lower = Kd_wall ! Store for next layer up.
      if (do_diag_Kd_BBL) Kd_BBL(i,j,K) = Kd_wall
    enddo ! k
  enddo ! i

end subroutine add_LOTW_BBL_diffusivity

!> This routine adds effects of mixed layer radiation to the layer diffusivities.
subroutine add_MLrad_diffusivity(h, fluxes, j, G, GV, US, CS, TKE_to_Kd, Kd_lay, Kd_int)
  type(ocean_grid_type),            intent(in)    :: G      !< The ocean's grid structure
  type(verticalGrid_type),          intent(in)    :: GV     !< The ocean's vertical grid structure
  type(unit_scale_type),            intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2]
  type(forcing),                    intent(in)    :: fluxes !< Surface fluxes structure
  integer,                          intent(in)    :: j      !< The j-index to work on
  type(set_diffusivity_CS),         pointer       :: CS     !< Diffusivity control structure
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: TKE_to_Kd !< The conversion rate between the TKE
                                                            !! TKE dissipated within  a layer and the
                                                            !! diapycnal diffusivity witin that layer,
                                                            !! usually (~Rho_0 / (G_Earth * dRho_lay))
                                                            !! [Z2 T-1 / Z3 T-3 = T2 Z-1 ~> s2 m-1]
  real, dimension(SZI_(G),SZK_(GV)), &
                          optional, intent(inout) :: Kd_lay !< The diapycnal diffusivity in layers [Z2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1), &
                          optional, intent(inout) :: Kd_int !< The diapycnal diffusivity at interfaces
                                                            !! [Z2 T-1 ~> m2 s-1].

! This routine adds effects of mixed layer radiation to the layer diffusivities.

  real, dimension(SZI_(G)) :: h_ml  ! Mixed layer thickness [Z ~> m].
  real, dimension(SZI_(G)) :: TKE_ml_flux ! Mixed layer TKE flux [Z3 T-3 ~> m3 s-3]
  real, dimension(SZI_(G)) :: I_decay ! A decay rate [Z-1 ~> m-1].
  real, dimension(SZI_(G)) :: Kd_mlr_ml ! Diffusivities associated with mixed layer radiation [Z2 T-1 ~> m2 s-1].

  real :: f_sq              ! The square of the local Coriolis parameter or a related variable [T-2 ~> s-2].
  real :: h_ml_sq           ! The square of the mixed layer thickness [Z2 ~> m2].
  real :: ustar_sq          ! ustar squared [Z2 T-2 ~> m2 s-2]
  real :: Kd_mlr            ! A diffusivity associated with mixed layer turbulence radiation [Z2 T-1 ~> m2 s-1].
  real :: C1_6              ! 1/6
  real :: Omega2            ! rotation rate squared [T-2 ~> s-2].
  real :: z1                ! layer thickness times I_decay [nondim]
  real :: dzL               ! thickness converted to heights [Z ~> m].
  real :: I_decay_len2_TKE  ! squared inverse decay lengthscale for
                            ! TKE, as used in the mixed layer code [Z-2 ~> m-2].
  real :: h_neglect         ! negligibly small thickness [Z ~> m].

  logical :: do_any, do_i(SZI_(G))
  integer :: i, k, is, ie, nz, kml
  is = G%isc ; ie = G%iec ; nz = GV%ke

  Omega2    = CS%omega**2
  C1_6      = 1.0 / 6.0
  kml       = GV%nkml
  h_neglect = GV%H_subroundoff*GV%H_to_Z

  if (.not.CS%ML_radiation) return

  do i=is,ie ; h_ml(i) = 0.0 ; do_i(i) = (G%mask2dT(i,j) > 0.5) ; enddo
  do k=1,kml ; do i=is,ie ; h_ml(i) = h_ml(i) + GV%H_to_Z*h(i,j,k) ; enddo ; enddo

  do i=is,ie ; if (do_i(i)) then
    if (CS%ML_omega_frac >= 1.0) then
      f_sq = 4.0 * Omega2
    else
      f_sq = 0.25 * ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                     (G%CoriolisBu(I,J-1)**2 + G%CoriolisBu(I-1,J)**2))
      if (CS%ML_omega_frac > 0.0) &
        f_sq = CS%ML_omega_frac * 4.0 * Omega2 + (1.0 - CS%ML_omega_frac) * f_sq
    endif

    ustar_sq = max(fluxes%ustar(i,j), CS%ustar_min)**2

    TKE_ml_flux(i) = (CS%mstar * CS%ML_rad_coeff) * (ustar_sq * (fluxes%ustar(i,j)))
    I_decay_len2_TKE = CS%TKE_decay**2 * (f_sq / ustar_sq)

    if (CS%ML_rad_TKE_decay) &
      TKE_ml_flux(i) = TKE_ml_flux(i) * exp(-h_ml(i) * sqrt(I_decay_len2_TKE))

    ! Calculate the inverse decay scale
    h_ml_sq = (CS%ML_rad_efold_coeff * (h_ml(i)+h_neglect))**2
    I_decay(i) = sqrt((I_decay_len2_TKE * h_ml_sq + 1.0) / h_ml_sq)

    ! Average the dissipation layer kml+1, using
    ! a more accurate Taylor series approximations for very thin layers.
    z1 = (GV%H_to_Z*h(i,j,kml+1)) * I_decay(i)
    if (z1 > 1e-5) then
      Kd_mlr = TKE_ml_flux(i) * TKE_to_Kd(i,kml+1) * (1.0 - exp(-z1))
    else
      Kd_mlr = TKE_ml_flux(i) * TKE_to_Kd(i,kml+1) * (z1 * (1.0 - z1 * (0.5 - C1_6 * z1)))
    endif
    Kd_mlr_ml(i) = min(Kd_mlr, CS%ML_rad_kd_max)
    TKE_ml_flux(i) = TKE_ml_flux(i) * exp(-z1)
  endif ; enddo

  if (present(Kd_lay)) then
    do k=1,kml+1 ; do i=is,ie ; if (do_i(i)) then
      Kd_lay(i,k) = Kd_lay(i,k) + Kd_mlr_ml(i)
    endif ; enddo ; enddo
  endif
  if (present(Kd_int)) then
    do K=2,kml+1 ; do i=is,ie ; if (do_i(i)) then
      Kd_int(i,K) = Kd_int(i,K) + Kd_mlr_ml(i)
    endif ; enddo ; enddo
    if (kml<=nz-1) then ; do i=is,ie ; if (do_i(i)) then
      Kd_int(i,Kml+2) = Kd_int(i,Kml+2) + 0.5 * Kd_mlr_ml(i)
    endif ; enddo ; endif
  endif

  do k=kml+2,nz-1
    do_any = .false.
    do i=is,ie ; if (do_i(i)) then
      dzL = GV%H_to_Z*h(i,j,k) ;  z1 = dzL*I_decay(i)
      if (CS%ML_Rad_bug) then
        ! These expresssions are dimensionally inconsistent. -RWH
        ! This is supposed to be the integrated energy deposited in the layer,
        ! not the average over the layer as in these expressions.
        if (z1 > 1e-5) then
          Kd_mlr = (TKE_ml_flux(i) * TKE_to_Kd(i,k)) * & ! Units of Z2 T-1
                   US%m_to_Z * ((1.0 - exp(-z1)) / dzL)  ! Units of m-1
        else
          Kd_mlr = (TKE_ml_flux(i) * TKE_to_Kd(i,k)) * &  ! Units of Z2 T-1
                   US%m_to_Z * (I_decay(i) * (1.0 - z1 * (0.5 - C1_6*z1))) ! Units of m-1
        endif
      else
        if (z1 > 1e-5) then
          Kd_mlr = (TKE_ml_flux(i) * TKE_to_Kd(i,k)) * (1.0 - exp(-z1))
        else
          Kd_mlr = (TKE_ml_flux(i) * TKE_to_Kd(i,k)) * (z1 * (1.0 - z1 * (0.5 - C1_6*z1)))
        endif
      endif
      Kd_mlr = min(Kd_mlr, CS%ML_rad_kd_max)
      if (present(Kd_lay)) then
        Kd_lay(i,k) = Kd_lay(i,k) + Kd_mlr
      endif
      if (present(Kd_int)) then
        Kd_int(i,K)   = Kd_int(i,K)   + 0.5 * Kd_mlr
        Kd_int(i,K+1) = Kd_int(i,K+1) + 0.5 * Kd_mlr
      endif

      TKE_ml_flux(i) = TKE_ml_flux(i) * exp(-z1)
      if (TKE_ml_flux(i) * I_decay(i) < 0.1 * CS%Kd_min * Omega2) then
        do_i(i) = .false.
      else ; do_any = .true. ; endif
    endif ; enddo
    if (.not.do_any) exit
  enddo

end subroutine add_MLrad_diffusivity

!> This subroutine calculates several properties related to bottom
!! boundary layer turbulence.
subroutine set_BBL_TKE(u, v, h, fluxes, visc, G, GV, US, CS, OBC)
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),  intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(forcing),            intent(in)    :: fluxes !< A structure of thermodynamic surface fluxes
  type(vertvisc_type),      intent(in)    :: visc !< Structure containing vertical viscosities, bottom
                                                  !! boundary layer properies, and related fields.
  type(set_diffusivity_CS), pointer       :: CS   !< Diffusivity control structure
  type(ocean_OBC_type), optional, pointer :: OBC  !< Open boundaries control structure.

  ! This subroutine calculates several properties related to bottom
  ! boundary layer turbulence.

  real, dimension(SZI_(G)) :: &
    htot          ! total thickness above or below a layer, or the
                  ! integrated thickness in the BBL [Z ~> m].

  real, dimension(SZIB_(G)) :: &
    uhtot, &      ! running integral of u in the BBL [Z L T-1 ~> m2 s-1]
    ustar, &      ! bottom boundary layer turbulence speed [Z T-1 ~> m s-1].
    u2_bbl        ! square of the mean zonal velocity in the BBL [L2 T-2 ~> m2 s-2]

  real :: vhtot(SZI_(G)) ! running integral of v in the BBL [Z L T-1 ~> m2 s-1]

  real, dimension(SZI_(G),SZJB_(G)) :: &
    vstar, & ! ustar at at v-points [Z T-1 ~> m s-1].
    v2_bbl   ! square of average meridional velocity in BBL [L2 T-2 ~> m2 s-2]

  real :: cdrag_sqrt  ! square root of the drag coefficient [nondim]
  real :: hvel        ! thickness at velocity points [Z ~> m].

  logical :: domore, do_i(SZI_(G))
  integer :: i, j, k, is, ie, js, je, nz
  integer :: l_seg
  logical :: local_open_u_BC, local_open_v_BC
  logical :: has_obc

  local_open_u_BC = .false.
  local_open_v_BC = .false.
  if (present(OBC)) then ; if (associated(OBC)) then
    local_open_u_BC = OBC%open_u_BCs_exist_globally
    local_open_v_BC = OBC%open_v_BCs_exist_globally
  endif ; endif

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) call MOM_error(FATAL,"set_BBL_TKE: "//&
         "Module must be initialized before it is used.")

  if (.not.CS%bottomdraglaw .or. (CS%BBL_effic<=0.0)) then
    if (associated(visc%ustar_BBL)) then
      do j=js,je ; do i=is,ie ; visc%ustar_BBL(i,j) = 0.0 ; enddo ; enddo
    endif
    if (associated(visc%TKE_BBL)) then
      do j=js,je ; do i=is,ie ; visc%TKE_BBL(i,j) = 0.0 ; enddo ; enddo
    endif
    return
  endif

  cdrag_sqrt = sqrt(CS%cdrag)

  !$OMP parallel default(shared) private(do_i,vhtot,htot,domore,hvel,uhtot,ustar,u2_bbl)
  !$OMP do
  do J=js-1,je
    ! Determine ustar and the square magnitude of the velocity in the
    ! bottom boundary layer. Together these give the TKE source and
    ! vertical decay scale.
    do i=is,ie ; if ((G%mask2dCv(i,J) > 0.5) .and. (cdrag_sqrt*visc%bbl_thick_v(i,J) > 0.0)) then
      do_i(i) = .true. ; vhtot(i) = 0.0 ; htot(i) = 0.0
      vstar(i,J) = visc%Kv_bbl_v(i,J) / (cdrag_sqrt*visc%bbl_thick_v(i,J))
    else
      do_i(i) = .false. ; vstar(i,J) = 0.0 ; htot(i) = 0.0
    endif ; enddo
    do k=nz,1,-1
      domore = .false.
      do i=is,ie ; if (do_i(i)) then
        ! Determine if grid point is an OBC
        has_obc = .false.
        if (local_open_v_BC) then
          l_seg = OBC%segnum_v(i,J)
          if (l_seg /= OBC_NONE) then
            has_obc = OBC%segment(l_seg)%open
          endif
        endif

        ! Compute h based on OBC state
        if (has_obc) then
          if (OBC%segment(l_seg)%direction == OBC_DIRECTION_N) then
            hvel = GV%H_to_Z*h(i,j,k)
          else
            hvel = GV%H_to_Z*h(i,j+1,k)
          endif
        else
          hvel = 0.5*GV%H_to_Z*(h(i,j,k) + h(i,j+1,k))
        endif

        if ((htot(i) + hvel) >= visc%bbl_thick_v(i,J)) then
          vhtot(i) = vhtot(i) + (visc%bbl_thick_v(i,J) - htot(i))*v(i,J,k)
          htot(i) = visc%bbl_thick_v(i,J)
          do_i(i) = .false.
        else
          vhtot(i) = vhtot(i) + hvel*v(i,J,k)
          htot(i) = htot(i) + hvel
          domore = .true.
        endif
      endif ; enddo
      if (.not.domore) exit
    enddo
    do i=is,ie ; if ((G%mask2dCv(i,J) > 0.5) .and. (htot(i) > 0.0)) then
      v2_bbl(i,J) = (vhtot(i)*vhtot(i))/(htot(i)*htot(i))
    else
      v2_bbl(i,J) = 0.0
    endif ; enddo
  enddo
  !$OMP do
  do j=js,je
    do I=is-1,ie ; if ((G%mask2dCu(I,j) > 0.5) .and. (cdrag_sqrt*visc%bbl_thick_u(I,j) > 0.0))  then
      do_i(I) = .true. ; uhtot(I) = 0.0 ; htot(I) = 0.0
      ustar(I) = visc%Kv_bbl_u(I,j) / (cdrag_sqrt*visc%bbl_thick_u(I,j))
    else
      do_i(I) = .false. ; ustar(I) = 0.0 ; htot(I) = 0.0
    endif ; enddo
    do k=nz,1,-1 ; domore = .false.
      do I=is-1,ie ; if (do_i(I)) then
        ! Determine if grid point is an OBC
        has_obc = .false.
        if (local_open_u_BC) then
          l_seg = OBC%segnum_u(I,j)
          if (l_seg /= OBC_NONE) then
            has_obc = OBC%segment(l_seg)%open
          endif
        endif

        ! Compute h based on OBC state
        if (has_obc) then
          if (OBC%segment(l_seg)%direction == OBC_DIRECTION_E) then
            hvel = GV%H_to_Z*h(i,j,k)
          else ! OBC_DIRECTION_W
            hvel = GV%H_to_Z*h(i+1,j,k)
          endif
        else
          hvel = 0.5*GV%H_to_Z*(h(i,j,k) + h(i+1,j,k))
        endif

        if ((htot(I) + hvel) >= visc%bbl_thick_u(I,j)) then
          uhtot(I) = uhtot(I) + (visc%bbl_thick_u(I,j) - htot(I))*u(I,j,k)
          htot(I) = visc%bbl_thick_u(I,j)
          do_i(I) = .false.
        else
          uhtot(I) = uhtot(I) + hvel*u(I,j,k)
          htot(I) = htot(I) + hvel
          domore = .true.
        endif
      endif ; enddo
      if (.not.domore) exit
    enddo
    do I=is-1,ie ; if ((G%mask2dCu(I,j) > 0.5) .and. (htot(i) > 0.0)) then
      u2_bbl(I) = (uhtot(I)*uhtot(I))/(htot(I)*htot(I))
    else
      u2_bbl(I) = 0.0
    endif ; enddo

    do i=is,ie
      visc%ustar_BBL(i,j) = sqrt(0.5*G%IareaT(i,j) * &
                ((G%areaCu(I-1,j)*(ustar(I-1)*ustar(I-1)) + &
                  G%areaCu(I,j)*(ustar(I)*ustar(I))) + &
                 (G%areaCv(i,J-1)*(vstar(i,J-1)*vstar(i,J-1)) + &
                  G%areaCv(i,J)*(vstar(i,J)*vstar(i,J))) ) )
      visc%TKE_BBL(i,j) = US%L_to_Z**2 * &
                 (((G%areaCu(I-1,j)*(ustar(I-1)*u2_bbl(I-1)) + &
                    G%areaCu(I,j) * (ustar(I)*u2_bbl(I))) + &
                   (G%areaCv(i,J-1)*(vstar(i,J-1)*v2_bbl(i,J-1)) + &
                    G%areaCv(i,J) * (vstar(i,J)*v2_bbl(i,J))) )*G%IareaT(i,j))
    enddo
  enddo
  !$OMP end parallel

end subroutine set_BBL_TKE

subroutine set_density_ratios(h, tv, kb, G, GV, US, CS, j, ds_dsp1, rho_0)
  type(ocean_grid_type),            intent(in)   :: G  !< The ocean's grid structure.
  type(verticalGrid_type),          intent(in)   :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)   :: h  !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),            intent(in)   :: tv !< Structure containing pointers to any
                                                       !! available thermodynamic fields; absent
                                                       !! fields have NULL ptrs.
  integer, dimension(SZI_(G)),      intent(in)   :: kb !< Index of lightest layer denser than the buffer
                                                       !! layer, or -1 without a bulk mixed layer.
  type(unit_scale_type),            intent(in)   :: US !< A dimensional unit scaling type
  type(set_diffusivity_CS),         pointer      :: CS !< Control structure returned by previous
                                                       !! call to diabatic_entrain_init.
  integer,                          intent(in)   :: j  !< Meridional index upon which to work.
  real, dimension(SZI_(G),SZK_(GV)), intent(out) :: ds_dsp1 !< Coordinate variable (sigma-2)
                                                       !! difference across an interface divided by
                                                       !! the difference across the interface below
                                                       !! it [nondim]
  real, dimension(SZI_(G),SZK_(GV)), &
                          optional, intent(in)   :: rho_0 !< Layer potential densities relative to
                                                       !! surface press [R ~> kg m-3].

  ! Local variables
  real :: g_R0                     ! g_R0 is a rescaled version of g/Rho [L2 Z-1 R-1 T-2 ~> m4 kg-1 s-2]
  real :: eps, tmp                 ! nondimensional temporary variables
  real :: a(SZK_(GV)), a_0(SZK_(GV)) ! nondimensional temporary variables
  real :: p_ref(SZI_(G))           ! an array of tv%P_Ref pressures [R L2 T-2 ~> Pa]
  real :: Rcv(SZI_(G),SZK_(GV))    ! coordinate density in the mixed and buffer layers [R ~> kg m-3]
  real :: I_Drho                   ! temporary variable [R-1 ~> m3 kg-1]

  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, k, k3, is, ie, nz, kmb
  is = G%isc ; ie = G%iec ; nz = GV%ke

  do k=2,nz-1
    if (GV%g_prime(k+1) /= 0.0) then
      do i=is,ie
        ds_dsp1(i,k) = GV%g_prime(k) / GV%g_prime(k+1)
      enddo
    else
      do i=is,ie
        ds_dsp1(i,k) = 1.
      enddo
    endif
  enddo

  if (CS%bulkmixedlayer) then
    g_R0 = GV%g_Earth / (GV%Rho0)
    kmb = GV%nk_rho_varies
    eps = 0.1
    do i=is,ie ; p_ref(i) = tv%P_Ref ; enddo
    EOSdom(:) = EOS_domain(G%HI)
    do k=1,kmb
      call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_ref, Rcv(:,k), tv%eqn_of_state, EOSdom)
    enddo
    do i=is,ie
      if (kb(i) <= nz-1) then
!   Set up appropriately limited ratios of the reduced gravities of the
! interfaces above and below the buffer layer and the next denser layer.
        k = kb(i)

        I_Drho = g_R0 / GV%g_prime(k+1)
        ! The indexing convention for a is appropriate for the interfaces.
        do k3=1,kmb
          a(k3+1) = (GV%Rlay(k) - Rcv(i,k3)) * I_Drho
        enddo
        if ((present(rho_0)) .and. (a(kmb+1) < 2.0*eps*ds_dsp1(i,k))) then
!   If the buffer layer nearly matches the density of the layer below in the
! coordinate variable (sigma-2), use the sigma-0-based density ratio if it is
! greater (and stable).
          if ((rho_0(i,k) > rho_0(i,kmb)) .and. &
              (rho_0(i,k+1) > rho_0(i,k))) then
            I_Drho = 1.0 / (rho_0(i,k+1)-rho_0(i,k))
            a_0(kmb+1) = min((rho_0(i,k)-rho_0(i,kmb)) * I_Drho, ds_dsp1(i,k))
            if (a_0(kmb+1) > a(kmb+1)) then
              do k3=2,kmb
                a_0(k3) = a_0(kmb+1) + (rho_0(i,kmb)-rho_0(i,k3-1)) * I_Drho
              enddo
              if (a(kmb+1) <= eps*ds_dsp1(i,k)) then
                do k3=2,kmb+1 ; a(k3) = a_0(k3) ; enddo
              else
! Alternative...  tmp = 0.5*(1.0 - cos(PI*(a(K2+1)/(eps*ds_dsp1(i,k)) - 1.0)) )
                tmp = a(kmb+1)/(eps*ds_dsp1(i,k)) - 1.0
                do k3=2,kmb+1 ; a(k3) = tmp*a(k3) + (1.0-tmp)*a_0(k3) ; enddo
              endif
            endif
          endif
        endif

        ds_dsp1(i,k) = MAX(a(kmb+1),1e-5)

        do k3=2,kmb
!           ds_dsp1(i,k3) = MAX(a(k3),1e-5)
          ! Deliberately treat convective instabilies of the upper mixed
          ! and buffer layers with respect to the deepest buffer layer as
          ! though they don't exist.  They will be eliminated by the upcoming
          ! call to the mixedlayer code anyway.
          ! The indexing convention is appropriate for the interfaces.
          ds_dsp1(i,k3) = MAX(a(k3),ds_dsp1(i,k))
        enddo
      endif ! (kb(i) <= nz-1)
    enddo ! I-loop.
  endif ! bulkmixedlayer

end subroutine set_density_ratios

subroutine set_diffusivity_init(Time, G, GV, US, param_file, diag, CS, int_tide_CSp, halo_TS, &
                                double_diffuse)
  type(time_type),          intent(in)    :: Time !< The current model time
  type(ocean_grid_type),    intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),    intent(in)    :: param_file !< A structure to parse for run-time
                                                  !! parameters.
  type(diag_ctrl), target,  intent(inout) :: diag !< A structure used to regulate diagnostic output.
  type(set_diffusivity_CS), pointer       :: CS   !< pointer set to point to the module control
                                                  !! structure.
  type(int_tide_CS),        pointer       :: int_tide_CSp !< A pointer to the internal tides control
                                                  !! structure
  integer,        optional, intent(out)   :: halo_TS !< The halo size of tracer points that must be
                                                  !! valid for the calculations in set_diffusivity.
  logical,        optional, intent(out)   :: double_diffuse !< If present, this indicates whether
                                                  !! some version of double diffusion is being used.

  ! Local variables
  real :: decay_length
  logical :: ML_use_omega
  logical :: default_2018_answers
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_set_diffusivity"  ! This module's name.
  real :: omega_frac_dflt
  logical :: Bryan_Lewis_diffusivity ! If true, the background diapycnal diffusivity uses
                                     ! the Bryan-Lewis (1979) style tanh profile.
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed

  if (associated(CS)) then
    call MOM_error(WARNING, "diabatic_entrain_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%diag => diag
  if (associated(int_tide_CSp))  CS%int_tide_CSp  => int_tide_CSp

  ! These default values always need to be set.
  CS%BBL_mixing_as_max = .true.
  CS%cdrag = 0.003 ; CS%BBL_effic = 0.0
  CS%bulkmixedlayer = (GV%nkml > 0)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, default=".")
  CS%inputdir = slasher(CS%inputdir)
  call get_param(param_file, mdl, "FLUX_RI_MAX", CS%FluxRi_max, &
                 "The flux Richardson number where the stratification is "//&
                 "large enough that N2 > omega2.  The full expression for "//&
                 "the Flux Richardson number is usually "//&
                 "FLUX_RI_MAX*N2/(N2+OMEGA2).", units="nondim", default=0.2)
  call get_param(param_file, mdl, "OMEGA", CS%omega, &
                 "The rotation rate of the earth.", units="s-1", default=7.2921e-5, scale=US%T_to_s)

  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  call get_param(param_file, mdl, "SET_DIFF_2018_ANSWERS", CS%answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)

  ! CS%use_tidal_mixing is set to True if an internal tidal dissipation scheme is to be used.
  CS%use_tidal_mixing = tidal_mixing_init(Time, G, GV, US, param_file, diag, &
                                          CS%tidal_mixing_CSp)

  call get_param(param_file, mdl, "ML_RADIATION", CS%ML_radiation, &
                 "If true, allow a fraction of TKE available from wind "//&
                 "work to penetrate below the base of the mixed layer "//&
                 "with a vertical decay scale determined by the minimum "//&
                 "of: (1) The depth of the mixed layer, (2) an Ekman "//&
                 "length scale.", default=.false.)
  if (CS%ML_radiation) then
    ! This give a minimum decay scale that is typically much less than Angstrom.
    CS%ustar_min = 2e-4 * CS%omega * (GV%Angstrom_Z + GV%H_subroundoff*GV%H_to_Z)

    call get_param(param_file, mdl, "ML_RAD_EFOLD_COEFF", CS%ML_rad_efold_coeff, &
                 "A coefficient that is used to scale the penetration "//&
                 "depth for turbulence below the base of the mixed layer. "//&
                 "This is only used if ML_RADIATION is true.", units="nondim", default=0.2)
    call get_param(param_file, mdl, "ML_RAD_BUG", CS%ML_rad_bug, &
                 "If true use code with a bug that reduces the energy available "//&
                 "in the transition layer by a factor of the inverse of the energy "//&
                 "deposition lenthscale (in m).", default=.false.)
    call get_param(param_file, mdl, "ML_RAD_KD_MAX", CS%ML_rad_kd_max, &
                 "The maximum diapycnal diffusivity due to turbulence "//&
                 "radiated from the base of the mixed layer. "//&
                 "This is only used if ML_RADIATION is true.", &
                 units="m2 s-1", default=1.0e-3, scale=US%m2_s_to_Z2_T)
    call get_param(param_file, mdl, "ML_RAD_COEFF", CS%ML_rad_coeff, &
                 "The coefficient which scales MSTAR*USTAR^3 to obtain "//&
                 "the energy available for mixing below the base of the "//&
                 "mixed layer. This is only used if ML_RADIATION is true.", &
                 units="nondim", default=0.2)
    call get_param(param_file, mdl, "ML_RAD_APPLY_TKE_DECAY", CS%ML_rad_TKE_decay, &
                 "If true, apply the same exponential decay to ML_rad as "//&
                 "is applied to the other surface sources of TKE in the "//&
                 "mixed layer code. This is only used if ML_RADIATION is true.", default=.true.)
    call get_param(param_file, mdl, "MSTAR", CS%mstar, &
                 "The ratio of the friction velocity cubed to the TKE "//&
                 "input to the mixed layer.", units="nondim", default=1.2)
    call get_param(param_file, mdl, "TKE_DECAY", CS%TKE_decay, &
                 "The ratio of the natural Ekman depth to the TKE decay scale.", &
                 units="nondim", default=2.5)
    call get_param(param_file, mdl, "ML_USE_OMEGA", ML_use_omega, &
                 "If true, use the absolute rotation rate instead of the "//&
                 "vertical component of rotation when setting the decay "//&
                 "scale for turbulence.", default=.false., do_not_log=.true.)
    omega_frac_dflt = 0.0
    if (ML_use_omega) then
      call MOM_error(WARNING, "ML_USE_OMEGA is depricated; use ML_OMEGA_FRAC=1.0 instead.")
      omega_frac_dflt = 1.0
    endif
    call get_param(param_file, mdl, "ML_OMEGA_FRAC", CS%ML_omega_frac, &
                   "When setting the decay scale for turbulence, use this "//&
                   "fraction of the absolute rotation rate blended with the "//&
                   "local value of f, as sqrt((1-of)*f^2 + of*4*omega^2).", &
                   units="nondim", default=omega_frac_dflt)
  endif

  call get_param(param_file, mdl, "BOTTOMDRAGLAW", CS%bottomdraglaw, &
                 "If true, the bottom stress is calculated with a drag "//&
                 "law of the form c_drag*|u|*u. The velocity magnitude "//&
                 "may be an assumed value or it may be based on the actual "//&
                 "velocity in the bottommost HBBL, depending on LINEAR_DRAG.", default=.true.)
  if  (CS%bottomdraglaw) then
    call get_param(param_file, mdl, "CDRAG", CS%cdrag, &
                 "The drag coefficient relating the magnitude of the "//&
                 "velocity field to the bottom stress. CDRAG is only used "//&
                 "if BOTTOMDRAGLAW is true.", units="nondim", default=0.003)
    call get_param(param_file, mdl, "BBL_EFFIC", CS%BBL_effic, &
                 "The efficiency with which the energy extracted by "//&
                 "bottom drag drives BBL diffusion.  This is only "//&
                 "used if BOTTOMDRAGLAW is true.", units="nondim", default=0.20)
    call get_param(param_file, mdl, "BBL_MIXING_MAX_DECAY", decay_length, &
                 "The maximum decay scale for the BBL diffusion, or 0 to allow the mixing "//&
                 "to penetrate as far as stratification and rotation permit.  The default "//&
                 "for now is 200 m. This is only used if BOTTOMDRAGLAW is true.", &
                 units="m", default=200.0, scale=US%m_to_Z)

    CS%IMax_decay = 0.0
    if (decay_length > 0.0) CS%IMax_decay = 1.0/decay_length
    call get_param(param_file, mdl, "BBL_MIXING_AS_MAX", CS%BBL_mixing_as_max, &
                 "If true, take the maximum of the diffusivity from the "//&
                 "BBL mixing and the other diffusivities. Otherwise, "//&
                 "diffusivity from the BBL_mixing is simply added.", &
                 default=.true.)
    call get_param(param_file, mdl, "USE_LOTW_BBL_DIFFUSIVITY", CS%use_LOTW_BBL_diffusivity, &
                 "If true, uses a simple, imprecise but non-coordinate dependent, model "//&
                 "of BBL mixing diffusivity based on Law of the Wall. Otherwise, uses "//&
                 "the original BBL scheme.", default=.false.)
    if (CS%use_LOTW_BBL_diffusivity) then
      call get_param(param_file, mdl, "LOTW_BBL_USE_OMEGA", CS%LOTW_BBL_use_omega, &
                 "If true, use the maximum of Omega and N for the TKE to diffusion "//&
                 "calculation. Otherwise, N is N.", default=.true.)
    endif
  else
    CS%use_LOTW_BBL_diffusivity = .false. ! This parameterization depends on a u* from viscous BBL
  endif
  CS%id_Kd_BBL = register_diag_field('ocean_model', 'Kd_BBL', diag%axesTi, Time, &
                 'Bottom Boundary Layer Diffusivity', 'm2 s-1', conversion=US%Z2_T_to_m2_s)
  call get_param(param_file, mdl, "SIMPLE_TKE_TO_KD", CS%simple_TKE_to_Kd, &
                 "If true, uses a simple estimate of Kd/TKE that will "//&
                 "work for arbitrary vertical coordinates. If false, "//&
                 "calculates Kd/TKE and bounds based on exact energetics "//&
                 "for an isopycnal layer-formulation.", default=.false.)

  ! set params related to the background mixing
  call bkgnd_mixing_init(Time, G, GV, US, param_file, CS%diag, CS%bkgnd_mixing_csp)

  call get_param(param_file, mdl, "KV", CS%Kv, &
                 "The background kinematic viscosity in the interior. "//&
                 "The molecular value, ~1e-6 m2 s-1, may be used.", &
                 units="m2 s-1", scale=US%m2_s_to_Z2_T, fail_if_missing=.true.)

  call get_param(param_file, mdl, "KD", CS%Kd, &
                 "The background diapycnal diffusivity of density in the "//&
                 "interior. Zero or the molecular value, ~1e-7 m2 s-1, "//&
                 "may be used.", default=0.0, units="m2 s-1", scale=US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "KD_MIN", CS%Kd_min, &
                 "The minimum diapycnal diffusivity.", &
                 units="m2 s-1", default=0.01*CS%Kd*US%Z2_T_to_m2_s, scale=US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "KD_MAX", CS%Kd_max, &
                 "The maximum permitted increment for the diapycnal "//&
                 "diffusivity from TKE-based parameterizations, or a negative "//&
                 "value for no limit.", units="m2 s-1", default=-1.0, scale=US%m2_s_to_Z2_T)
  if (CS%simple_TKE_to_Kd .and. CS%Kd_max<=0.) call MOM_error(FATAL, &
         "set_diffusivity_init: To use SIMPLE_TKE_TO_KD, KD_MAX must be set to >0.")
  call get_param(param_file, mdl, "KD_ADD", CS%Kd_add, &
                 "A uniform diapycnal diffusivity that is added "//&
                 "everywhere without any filtering or scaling.", &
                 units="m2 s-1", default=0.0, scale=US%m2_s_to_Z2_T)
  if (CS%use_LOTW_BBL_diffusivity .and. CS%Kd_max<=0.) call MOM_error(FATAL, &
                 "set_diffusivity_init: KD_MAX must be set (positive) when "// &
                 "USE_LOTW_BBL_DIFFUSIVITY=True.")
  call get_param(param_file, mdl, "KD_SMOOTH", CS%Kd_smooth, &
                 "A diapycnal diffusivity that is used to interpolate "//&
                 "more sensible values of T & S into thin layers.", &
                 units="m2 s-1", default=1.0e-6, scale=US%m2_s_to_Z2_T)

  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)

  call get_param(param_file, mdl, "USER_CHANGE_DIFFUSIVITY", CS%user_change_diff, &
                 "If true, call user-defined code to change the diffusivity.", default=.false.)

  call get_param(param_file, mdl, "DISSIPATION_MIN", CS%dissip_min, &
                 "The minimum dissipation by which to determine a lower "//&
                 "bound of Kd (a floor).", &
                 units="W m-3", default=0.0, scale=US%W_m2_to_RZ3_T3*US%Z_to_m)
  call get_param(param_file, mdl, "DISSIPATION_N0", CS%dissip_N0, &
                 "The intercept when N=0 of the N-dependent expression "//&
                 "used to set a minimum dissipation by which to determine "//&
                 "a lower bound of Kd (a floor): A in eps_min = A + B*N.", &
                 units="W m-3", default=0.0, scale=US%W_m2_to_RZ3_T3*US%Z_to_m)
  call get_param(param_file, mdl, "DISSIPATION_N1", CS%dissip_N1, &
                 "The coefficient multiplying N, following Gargett, used to "//&
                 "set a minimum dissipation by which to determine a lower "//&
                 "bound of Kd (a floor): B in eps_min = A + B*N", &
                 units="J m-3", default=0.0, scale=US%W_m2_to_RZ3_T3*US%Z_to_m*US%s_to_T)
  call get_param(param_file, mdl, "DISSIPATION_KD_MIN", CS%dissip_Kd_min, &
                 "The minimum vertical diffusivity applied as a floor.", &
                 units="m2 s-1", default=0.0, scale=US%m2_s_to_Z2_T)

  CS%limit_dissipation = (CS%dissip_min>0.) .or. (CS%dissip_N1>0.) .or. &
                         (CS%dissip_N0>0.) .or. (CS%dissip_Kd_min>0.)
  CS%dissip_N2 = 0.0
  if (CS%FluxRi_max > 0.0) &
    CS%dissip_N2 = CS%dissip_Kd_min * GV%Rho0 / CS%FluxRi_max

  CS%id_Kd_bkgnd = register_diag_field('ocean_model', 'Kd_bkgnd', diag%axesTi, Time, &
      'Background diffusivity added by MOM_bkgnd_mixing module', 'm2/s', conversion=US%Z2_T_to_m2_s)
  CS%id_Kv_bkgnd = register_diag_field('ocean_model', 'Kv_bkgnd', diag%axesTi, Time, &
      'Background viscosity added by MOM_bkgnd_mixing module', 'm2/s', conversion=US%Z2_T_to_m2_s)

  CS%id_Kd_layer = register_diag_field('ocean_model', 'Kd_layer', diag%axesTL, Time, &
      'Diapycnal diffusivity of layers (as set)', 'm2 s-1', conversion=US%Z2_T_to_m2_s)

  if (CS%use_tidal_mixing) then
    CS%id_Kd_Work = register_diag_field('ocean_model', 'Kd_Work', diag%axesTL, Time, &
         'Work done by Diapycnal Mixing', 'W m-2', conversion=US%RZ3_T3_to_W_m2)
    CS%id_maxTKE = register_diag_field('ocean_model', 'maxTKE', diag%axesTL, Time, &
           'Maximum layer TKE', 'm3 s-3', conversion=(US%Z_to_m**3*US%s_to_T**3))
    CS%id_TKE_to_Kd = register_diag_field('ocean_model', 'TKE_to_Kd', diag%axesTL, Time, &
           'Convert TKE to Kd', 's2 m', conversion=US%Z2_T_to_m2_s*(US%m_to_Z**3*US%T_to_s**3))
    CS%id_N2 = register_diag_field('ocean_model', 'N2', diag%axesTi, Time, &
         'Buoyancy frequency squared', 's-2', conversion=US%s_to_T**2, cmor_field_name='obvfsq', &
          cmor_long_name='Square of seawater buoyancy frequency', &
          cmor_standard_name='square_of_brunt_vaisala_frequency_in_sea_water')
  endif

  if (CS%user_change_diff) &
    CS%id_Kd_user = register_diag_field('ocean_model', 'Kd_user', diag%axesTi, Time, &
         'User-specified Extra Diffusivity', 'm2 s-1', conversion=US%Z2_T_to_m2_s)

  call get_param(param_file, mdl, "DOUBLE_DIFFUSION", CS%double_diffusion, &
                 "If true, increase diffusivites for temperature or salinity based on the "//&
                 "double-diffusive parameterization described in Large et al. (1994).", &
                 default=.false.)

  if (CS%double_diffusion) then
    call get_param(param_file, mdl, "MAX_RRHO_SALT_FINGERS", CS%Max_Rrho_salt_fingers, &
                 "Maximum density ratio for salt fingering regime.", &
                 default=2.55, units="nondim")
    call get_param(param_file, mdl, "MAX_SALT_DIFF_SALT_FINGERS", CS%Max_salt_diff_salt_fingers, &
                 "Maximum salt diffusivity for salt fingering regime.", &
                 default=1.e-4, units="m2 s-1", scale=US%m2_s_to_Z2_T)
    call get_param(param_file, mdl, "KV_MOLECULAR", CS%Kv_molecular, &
                 "Molecular viscosity for calculation of fluxes under double-diffusive "//&
                 "convection.", default=1.5e-6, units="m2 s-1", scale=US%m2_s_to_Z2_T)
    ! The default molecular viscosity follows the CCSM4.0 and MOM4p1 defaults.
  endif ! old double-diffusion

  if (CS%user_change_diff) then
    call user_change_diff_init(Time, G, GV, US, param_file, diag, CS%user_change_diff_CSp)
  endif

  call get_param(param_file, mdl, "BRYAN_LEWIS_DIFFUSIVITY", Bryan_Lewis_diffusivity, &
                 "If true, use a Bryan & Lewis (JGR 1979) like tanh "//&
                 "profile of background diapycnal diffusivity with depth. "//&
                 "This is done via CVMix.", default=.false., do_not_log=.true.)
  if (CS%use_tidal_mixing .and. Bryan_Lewis_diffusivity) &
    call MOM_error(FATAL,"MOM_Set_Diffusivity: "// &
         "Bryan-Lewis and internal tidal dissipation are both enabled. Choose one.")

  CS%useKappaShear = kappa_shear_init(Time, G, GV, US, param_file, CS%diag, CS%kappaShear_CSp)
  CS%Vertex_Shear = kappa_shear_at_vertex(param_file)

  if (CS%useKappaShear) &
    id_clock_kappaShear = cpu_clock_id('(Ocean kappa_shear)', grain=CLOCK_MODULE)

  ! CVMix shear-driven mixing
  CS%use_CVMix_shear = CVMix_shear_init(Time, G, GV, US, param_file, CS%diag, CS%CVMix_shear_csp)

  ! CVMix double diffusion mixing
  CS%use_CVMix_ddiff = CVMix_ddiff_init(Time, G, GV, US, param_file, CS%diag, CS%CVMix_ddiff_csp)
  if (CS%use_CVMix_ddiff) &
    id_clock_CVMix_ddiff = cpu_clock_id('(Double diffusion via CVMix)', grain=CLOCK_MODULE)

  if (CS%double_diffusion .and. CS%use_CVMix_ddiff) then
    call MOM_error(FATAL, 'set_diffusivity_init: '// &
           'Multiple double-diffusion options selected (DOUBLE_DIFFUSION and'//&
           'USE_CVMIX_DDIFF), please disable all but one option to proceed.')
  endif

  if (CS%double_diffusion .or. CS%use_CVMix_ddiff) then
    CS%id_KT_extra = register_diag_field('ocean_model', 'KT_extra', diag%axesTi, Time, &
         'Double-diffusive diffusivity for temperature', 'm2 s-1', conversion=US%Z2_T_to_m2_s)
    CS%id_KS_extra = register_diag_field('ocean_model', 'KS_extra', diag%axesTi, Time, &
         'Double-diffusive diffusivity for salinity', 'm2 s-1', conversion=US%Z2_T_to_m2_s)
  endif
  if (CS%use_CVMix_ddiff) then
    CS%id_R_rho = register_diag_field('ocean_model', 'R_rho', diag%axesTi, Time, &
         'Double-diffusion density ratio', 'nondim')
  endif

  if (present(halo_TS)) then
    halo_TS = 0
    if (CS%Vertex_Shear) halo_TS = 1
  endif

  if (present(double_diffuse)) then
    double_diffuse = (CS%double_diffusion .or. CS%use_CVMix_ddiff)
  endif

end subroutine set_diffusivity_init

!> Clear pointers and dealocate memory
subroutine set_diffusivity_end(CS)
  type(set_diffusivity_CS), pointer :: CS !< Control structure for this module

  if (.not.associated(CS)) return

  call bkgnd_mixing_end(CS%bkgnd_mixing_csp)

  if (CS%use_tidal_mixing) call tidal_mixing_end(CS%tidal_mixing_CSp)

  if (CS%user_change_diff) call user_change_diff_end(CS%user_change_diff_CSp)

  if (CS%use_CVMix_shear)  call CVMix_shear_end(CS%CVMix_shear_csp)

  if (CS%use_CVMix_ddiff)  call CVMix_ddiff_end(CS%CVMix_ddiff_csp)

  if (associated(CS)) deallocate(CS)

end subroutine set_diffusivity_end

end module MOM_set_diffusivity
