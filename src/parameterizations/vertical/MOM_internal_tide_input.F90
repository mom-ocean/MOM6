!> Calculates energy input to the internal tides
module MOM_int_tide_input

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,        only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,    only : diag_ctrl, query_averaging_enabled
use MOM_diag_mediator,    only : disable_averaging, enable_averages
use MOM_diag_mediator,    only : safe_alloc_ptr, post_data, register_diag_field
use MOM_debugging,        only : hchksum
use MOM_error_handler,    only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,      only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,     only : forcing
use MOM_grid,             only : ocean_grid_type
use MOM_io,               only : slasher, vardesc, MOM_read_data
use MOM_isopycnal_slopes, only : vert_fill_TS
use MOM_time_manager,     only : time_type, set_time, operator(+), operator(<=)
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : thermo_var_ptrs, vertvisc_type, p3d
use MOM_verticalGrid,     only : verticalGrid_type
use MOM_EOS,              only : calculate_density, calculate_density_derivs, EOS_domain

implicit none ; private

#include <MOM_memory.h>

public set_int_tide_input, int_tide_input_init, int_tide_input_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> This control structure holds parameters that regulate internal tide energy inputs.
type, public :: int_tide_input_CS ; private
  logical :: debug      !< If true, write verbose checksums for debugging.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                        !! regulate the timing of diagnostic output.
  real :: TKE_itide_max !< Maximum Internal tide conversion
                        !! available to mix above the BBL [R Z3 T-3 ~> W m-2]
  real :: kappa_fill    !< Vertical diffusivity used to interpolate sensible values
                        !! of T & S into thin layers [Z2 T-1 ~> m2 s-1].

  real, allocatable, dimension(:,:) :: TKE_itidal_coef
            !< The time-invariant field that enters the TKE_itidal input calculation [R Z3 T-2 ~> J m-2].
  character(len=200) :: inputdir !< The directory for input files.

  logical :: int_tide_source_test    !< If true, apply an arbitrary generation site
                                     !! for internal tide testing (BDM)
  type(time_type) :: time_max_source !< A time for use in testing internal tides
  real    :: int_tide_source_x       !< X Location of generation site
                                     !! for internal tide for testing (BDM)
  real    :: int_tide_source_y       !< Y Location of generation site
                                     !! for internal tide for testing (BDM)


  !>@{ Diagnostic IDs
  integer :: id_TKE_itidal_itide = -1, id_Nb = -1, id_N2_bot = -1
  !>@}
end type int_tide_input_CS

!> This type is used to exchange fields related to the internal tides.
type, public :: int_tide_input_type
  real, allocatable, dimension(:,:) :: &
    TKE_itidal_input, & !< The internal tide TKE input at the bottom of the ocean [R Z3 T-3 ~> W m-2].
    h2, &               !< The squared topographic roughness height [Z2 ~> m2].
    tideamp, &          !< The amplitude of the tidal velocities [Z T-1 ~> m s-1].
    Nb                  !< The bottom stratification [T-1 ~> s-1].
end type int_tide_input_type

contains

!> Sets the model-state dependent internal tide energy sources.
subroutine set_int_tide_input(u, v, h, tv, fluxes, itide, dt, G, GV, US, CS)
  type(ocean_grid_type),                     intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< The ocean's vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u  !< The zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v  !< The meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(in)    :: tv !< A structure containing pointers to the
                                                                 !! thermodynamic fields
  type(forcing),                             intent(in)    :: fluxes !< A structure of thermodynamic surface fluxes
  type(int_tide_input_type),                 intent(inout) :: itide !< A structure containing fields related
                                                                 !! to the internal tide sources.
  real,                                      intent(in)    :: dt !< The time increment [T ~> s].
  type(int_tide_input_CS),                   pointer       :: CS !< This module's control structure.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    N2_bot        ! The bottom squared buoyancy frequency [T-2 ~> s-2].

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    T_f, S_f      ! The temperature and salinity in [degC] and [ppt] with the values in
                  ! the massless layers filled vertically by diffusion.
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  logical :: avg_enabled  ! for testing internal tides (BDM)
  type(time_type) :: time_end        !< For use in testing internal tides (BDM)

  integer :: i, j, k, is, ie, js, je, nz, isd, ied, jsd, jed

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) call MOM_error(FATAL,"set_diffusivity: "//&
         "Module must be initialized before it is used.")

  use_EOS = associated(tv%eqn_of_state)

  ! Smooth the properties through massless layers.
  if (use_EOS) then
    call vert_fill_TS(h, tv%T, tv%S, CS%kappa_fill*dt, T_f, S_f, G, GV, larger_h_denom=.true.)
  endif

  call find_N2_bottom(h, tv, T_f, S_f, itide%h2, fluxes, G, GV, US, N2_bot)

  avg_enabled = query_averaging_enabled(CS%diag, time_end=time_end)

  !$OMP parallel do default(shared)
  do j=js,je ; do i=is,ie
    itide%Nb(i,j) = G%mask2dT(i,j) * sqrt(N2_bot(i,j))
    itide%TKE_itidal_input(i,j) = min(CS%TKE_itidal_coef(i,j)*itide%Nb(i,j), CS%TKE_itide_max)
  enddo ; enddo

  if (CS%int_tide_source_test) then
    itide%TKE_itidal_input(:,:) = 0.0
    if (time_end <= CS%time_max_source) then
      do j=js,je ; do i=is,ie
        ! Input  an arbitrary energy point source.id_
        if (((G%geoLonCu(I-1,j)-CS%int_tide_source_x) * (G%geoLonBu(I,j)-CS%int_tide_source_x) <= 0.0) .and. &
            ((G%geoLatCv(i,J-1)-CS%int_tide_source_y) * (G%geoLatCv(i,j)-CS%int_tide_source_y) <= 0.0)) then
          itide%TKE_itidal_input(i,j) = 1.0*US%kg_m3_to_R*US%m_to_Z**3*US%T_to_s**3
        endif
      enddo ; enddo
    endif
  endif

  if (CS%debug) then
    call hchksum(N2_bot,"N2_bot",G%HI,haloshift=0, scale=US%s_to_T**2)
    call hchksum(itide%TKE_itidal_input,"TKE_itidal_input",G%HI,haloshift=0, &
                 scale=US%RZ3_T3_to_W_m2)
  endif

  call enable_averages(dt, time_end, CS%diag)

  if (CS%id_TKE_itidal_itide > 0) call post_data(CS%id_TKE_itidal_itide, itide%TKE_itidal_input, CS%diag)
  if (CS%id_Nb > 0) call post_data(CS%id_Nb, itide%Nb, CS%diag)
  if (CS%id_N2_bot > 0 ) call post_data(CS%id_N2_bot, N2_bot, CS%diag)

  call disable_averaging(CS%diag)

end subroutine set_int_tide_input

!> Estimates the near-bottom buoyancy frequency (N^2).
subroutine find_N2_bottom(h, tv, T_f, S_f, h2, fluxes, G, GV, US, N2_bot)
  type(ocean_grid_type),                    intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),                  intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                    intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                    intent(in)  :: tv   !< A structure containing pointers to the
                                                                !! thermodynamic fields
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)  :: T_f  !< Temperature after vertical filtering to
                                                                !! smooth out the values in thin layers [degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)  :: S_f  !< Salinity after vertical filtering to
                                                                !! smooth out the values in thin layers [ppt].
  real, dimension(SZI_(G),SZJ_(G)),         intent(in)  :: h2   !< Bottom topographic roughness [Z2 ~> m2].
  type(forcing),                            intent(in)  :: fluxes !< A structure of thermodynamic surface fluxes
  type(int_tide_input_CS),                  pointer     :: CS    !<  This module's control structure.
  real, dimension(SZI_(G),SZJ_(G)),         intent(out) :: N2_bot !< The squared buoyancy freqency at the
                                                                 !! ocean bottom [T-2 ~> s-2].
  ! Local variables
  real, dimension(SZI_(G),SZK_(G)+1) :: &
    dRho_int      ! The unfiltered density differences across interfaces [R ~> kg m-3].
  real, dimension(SZI_(G)) :: &
    pres, &       ! The pressure at each interface [R L2 T-2 ~> Pa].
    Temp_int, &   ! The temperature at each interface [degC].
    Salin_int, &  ! The salinity at each interface [ppt].
    drho_bot, &   ! The density difference at the bottom of a layer [R ~> kg m-3]
    h_amp, &      ! The amplitude of topographic roughness [Z ~> m].
    hb, &         ! The depth below a layer [Z ~> m].
    z_from_bot, & ! The height of a layer center above the bottom [Z ~> m].
    dRho_dT, &    ! The partial derivative of density with temperature [R degC-1 ~> kg m-3 degC-1]
    dRho_dS       ! The partial derivative of density with salinity [R ppt-1 ~> kg m-3 ppt-1].

  real :: dz_int  ! The thickness associated with an interface [Z ~> m].
  real :: G_Rho0  ! The gravitation acceleration divided by the Boussinesq
                  ! density [Z T-2 R-1 ~> m4 s-2 kg-1].
  logical :: do_i(SZI_(G)), do_any
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  G_Rho0 = (US%L_to_Z**2*GV%g_Earth) / GV%Rho0
  EOSdom(:) = EOS_domain(G%HI)

  ! Find the (limited) density jump across each interface.
  do i=is,ie
    dRho_int(i,1) = 0.0 ; dRho_int(i,nz+1) = 0.0
  enddo
!$OMP parallel do default(none) shared(is,ie,js,je,nz,tv,fluxes,G,GV,US,h,T_f,S_f, &
!$OMP                                  h2,N2_bot,G_Rho0,EOSdom) &
!$OMP                          private(pres,Temp_Int,Salin_Int,dRho_dT,dRho_dS, &
!$OMP                                  hb,dRho_bot,z_from_bot,do_i,h_amp,       &
!$OMP                                  do_any,dz_int) &
!$OMP                     firstprivate(dRho_int)
  do j=js,je
    if (associated(tv%eqn_of_state)) then
      if (associated(fluxes%p_surf)) then
        do i=is,ie ; pres(i) = fluxes%p_surf(i,j) ; enddo
      else
        do i=is,ie ; pres(i) = 0.0 ; enddo
      endif
      do K=2,nz
        do i=is,ie
          pres(i) = pres(i) + (GV%g_Earth*GV%H_to_RZ)*h(i,j,k-1)
          Temp_Int(i) = 0.5 * (T_f(i,j,k) + T_f(i,j,k-1))
          Salin_Int(i) = 0.5 * (S_f(i,j,k) + S_f(i,j,k-1))
        enddo
        call calculate_density_derivs(Temp_int, Salin_int, pres, dRho_dT(:), dRho_dS(:), &
                                      tv%eqn_of_state, EOSdom)
        do i=is,ie
          dRho_int(i,K) = max(dRho_dT(i)*(T_f(i,j,k) - T_f(i,j,k-1)) + &
                              dRho_dS(i)*(S_f(i,j,k) - S_f(i,j,k-1)), 0.0)
        enddo
      enddo
    else
      do K=2,nz ; do i=is,ie
        dRho_int(i,K) = (GV%Rlay(k) - GV%Rlay(k-1))
      enddo ; enddo
    endif

    ! Find the bottom boundary layer stratification.
    do i=is,ie
      hb(i) = 0.0 ; dRho_bot(i) = 0.0
      z_from_bot(i) = 0.5*GV%H_to_Z*h(i,j,nz)
      do_i(i) = (G%mask2dT(i,j) > 0.5)
      h_amp(i) = sqrt(h2(i,j))
    enddo

    do k=nz,2,-1
      do_any = .false.
      do i=is,ie ; if (do_i(i)) then
        dz_int = 0.5*GV%H_to_Z*(h(i,j,k) + h(i,j,k-1))
        z_from_bot(i) = z_from_bot(i) + dz_int ! middle of the layer above

        hb(i) = hb(i) + dz_int
        dRho_bot(i) = dRho_bot(i) + dRho_int(i,K)

        if (z_from_bot(i) > h_amp(i)) then
          if (k>2) then
            ! Always include at least one full layer.
            hb(i) = hb(i) + 0.5*GV%H_to_Z*(h(i,j,k-1) + h(i,j,k-2))
            dRho_bot(i) = dRho_bot(i) + dRho_int(i,K-1)
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
        N2_bot(i,j) = (G_Rho0 * dRho_bot(i)) / hb(i)
      else ;  N2_bot(i,j) = 0.0 ; endif
    enddo
  enddo

end subroutine find_N2_bottom

!> Initializes the data related to the internal tide input module
subroutine int_tide_input_init(Time, G, GV, US, param_file, diag, CS, itide)
  type(time_type),           intent(in)    :: Time !< The current model time
  type(ocean_grid_type),     intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),     intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),     intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl),   target, intent(inout) :: diag !< structure used to regulate diagnostic output.
  type(int_tide_input_CS),   pointer       :: CS   !< This module's control structure, which is initialized here.
  type(int_tide_input_type), pointer       :: itide !< A structure containing fields related
                                                   !! to the internal tide sources.
  ! Local variables
  type(vardesc) :: vd
  logical :: read_tideamp
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_int_tide_input"  ! This module's name.
  character(len=20)  :: tmpstr
  character(len=200) :: filename, tideamp_file, h2_file

  real :: mask_itidal        ! A multiplicative land mask, 0 or 1 [nondim]
  real :: max_frac_rough     ! The fraction relating the maximum topographic roughness
                             ! to the mean depth [nondim]
  real :: utide              ! constant tidal amplitude [L T-1 ~> m s-1] to be used if
                             ! tidal amplitude file is not present.
  real :: kappa_h2_factor    ! factor for the product of wavenumber * rms sgs height [nondim].
  real :: kappa_itides       ! topographic wavenumber and non-dimensional scaling [L-1 ~> m-1]
  real :: min_zbot_itides    ! Minimum ocean depth for internal tide conversion [Z ~> m].
  integer :: tlen_days       !< Time interval from start for adding wave source
                             !! for testing internal tides (BDM)
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed

  if (associated(CS)) then
    call MOM_error(WARNING, "int_tide_input_init called with an associated "// &
                            "control structure.")
    return
  endif
  if (associated(itide)) then
    call MOM_error(WARNING, "int_tide_input_init called with an associated "// &
                            "internal tide input type.")
    return
  endif
  allocate(CS)
  allocate(itide)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, default=".")
  CS%inputdir = slasher(CS%inputdir)

  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false., do_not_log=.true.)

  call get_param(param_file, mdl, "MIN_ZBOT_ITIDES", min_zbot_itides, &
               "Turn off internal tidal dissipation when the total "//&
               "ocean depth is less than this value.", units="m", default=0.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "KD_SMOOTH", CS%kappa_fill, &
                 "A diapycnal diffusivity that is used to interpolate "//&
                 "more sensible values of T & S into thin layers.", &
                 units="m2 s-1", default=1.0e-6, scale=US%m2_s_to_Z2_T)

  call get_param(param_file, mdl, "UTIDE", utide, &
               "The constant tidal amplitude used with INT_TIDE_DISSIPATION.", &
               units="m s-1", default=0.0, scale=US%m_s_to_L_T)

  allocate(itide%Nb(isd:ied,jsd:jed))  ; itide%Nb(:,:) = 0.0
  allocate(itide%h2(isd:ied,jsd:jed))  ; itide%h2(:,:) = 0.0
  allocate(itide%TKE_itidal_input(isd:ied,jsd:jed)) ; itide%TKE_itidal_input(:,:) = 0.0
  allocate(itide%tideamp(isd:ied,jsd:jed)) ; itide%tideamp(:,:) = utide
  allocate(CS%TKE_itidal_coef(isd:ied,jsd:jed)) ; CS%TKE_itidal_coef(:,:) = 0.0

  call get_param(param_file, mdl, "KAPPA_ITIDES", kappa_itides, &
               "A topographic wavenumber used with INT_TIDE_DISSIPATION. "//&
               "The default is 2pi/10 km, as in St.Laurent et al. 2002.", &
               units="m-1", default=8.e-4*atan(1.0), scale=US%L_to_m)

  call get_param(param_file, mdl, "KAPPA_H2_FACTOR", kappa_h2_factor, &
               "A scaling factor for the roughness amplitude with n"//&
               "INT_TIDE_DISSIPATION.",  units="nondim", default=1.0)
  call get_param(param_file, mdl, "TKE_ITIDE_MAX", CS%TKE_itide_max, &
               "The maximum internal tide energy source available to mix "//&
               "above the bottom boundary layer with INT_TIDE_DISSIPATION.", &
               units="W m-2", default=1.0e3, scale=US%W_m2_to_RZ3_T3)

  call get_param(param_file, mdl, "READ_TIDEAMP", read_tideamp, &
               "If true, read a file (given by TIDEAMP_FILE) containing "//&
               "the tidal amplitude with INT_TIDE_DISSIPATION.", default=.false.)
  if (read_tideamp) then
    call get_param(param_file, mdl, "TIDEAMP_FILE", tideamp_file, &
               "The path to the file containing the spatially varying "//&
               "tidal amplitudes with INT_TIDE_DISSIPATION.", default="tideamp.nc")
    filename = trim(CS%inputdir) // trim(tideamp_file)
    call log_param(param_file, mdl, "INPUTDIR/TIDEAMP_FILE", filename)
    call MOM_read_data(filename, 'tideamp', itide%tideamp, G%domain, timelevel=1, scale=US%m_s_to_L_T)
  endif

  call get_param(param_file, mdl, "H2_FILE", h2_file, &
               "The path to the file containing the sub-grid-scale "//&
               "topographic roughness amplitude with INT_TIDE_DISSIPATION.", &
               fail_if_missing=.true.)
  filename = trim(CS%inputdir) // trim(h2_file)
  call log_param(param_file, mdl, "INPUTDIR/H2_FILE", filename)
  call MOM_read_data(filename, 'h2', itide%h2, G%domain, timelevel=1, scale=US%m_to_Z**2)

  call get_param(param_file, mdl, "FRACTIONAL_ROUGHNESS_MAX", max_frac_rough, &
                 "The maximum topographic roughness amplitude as a fraction of the mean depth, "//&
                 "or a negative value for no limitations on roughness.", &
                 units="nondim", default=0.1)

  ! The following parameters are used in testing the internal tide code.
  call get_param(param_file, mdl, "INTERNAL_TIDE_SOURCE_TEST", CS%int_tide_source_test, &
                 "If true, apply an arbitrary generation site for internal tide testing", &
                 default=.false.)
  if (CS%int_tide_source_test)then
    call get_param(param_file, mdl, "INTERNAL_TIDE_SOURCE_X", CS%int_tide_source_x, &
                 "X Location of generation site for internal tide", default=1.)
    call get_param(param_file, mdl, "INTERNAL_TIDE_SOURCE_Y", CS%int_tide_source_y, &
                 "Y Location of generation site for internal tide", default=1.)
    call get_param(param_file, mdl, "INTERNAL_TIDE_SOURCE_TLEN_DAYS", tlen_days, &
                 "Time interval from start of experiment for adding wave source", &
                 units="days", default=0)
    CS%time_max_source = Time + set_time(0, days=tlen_days)
  endif

  do j=js,je ; do i=is,ie
    mask_itidal = 1.0
    if (G%bathyT(i,j) < min_zbot_itides) mask_itidal = 0.0

    itide%tideamp(i,j) = itide%tideamp(i,j) * mask_itidal * G%mask2dT(i,j)

    ! Restrict rms topo to a fraction (often 10 percent) of the column depth.
    if (max_frac_rough >= 0.0) &
      itide%h2(i,j) = min((max_frac_rough*G%bathyT(i,j))**2, itide%h2(i,j))

    ! Compute the fixed part of internal tidal forcing; units are [R Z3 T-2 ~> J m-2] here.
    CS%TKE_itidal_coef(i,j) = 0.5*US%L_to_Z*kappa_h2_factor*GV%Rho0*&
         kappa_itides * itide%h2(i,j) * itide%tideamp(i,j)**2
  enddo ; enddo


  CS%id_TKE_itidal_itide = register_diag_field('ocean_model','TKE_itidal_itide',diag%axesT1,Time, &
      'Internal Tide Driven Turbulent Kinetic Energy', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2)

  CS%id_Nb = register_diag_field('ocean_model','Nb_itide',diag%axesT1,Time, &
       'Bottom Buoyancy Frequency', 's-1', conversion=US%s_to_T)

  CS%id_N2_bot = register_diag_field('ocean_model','N2_b_itide',diag%axesT1,Time, &
       'Bottom Buoyancy frequency squared', 's-2', conversion=US%s_to_T**2)

end subroutine int_tide_input_init

!> Deallocates any memory related to the internal tide input module.
subroutine int_tide_input_end(CS)
  type(int_tide_input_CS), pointer :: CS !< This module's control structure, which is deallocated here.

  if (associated(CS)) deallocate(CS)

end subroutine int_tide_input_end

end module MOM_int_tide_input
