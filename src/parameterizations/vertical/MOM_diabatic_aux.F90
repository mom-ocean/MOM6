!> Provides functions for some diabatic processes such as fraxil, brine rejection,
!! tendency due to surface flux divergence.
module MOM_diabatic_aux

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,     only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_EOS,           only : calculate_density, calculate_TFreeze, EOS_domain
use MOM_EOS,           only : calculate_specific_vol_derivs, calculate_density_derivs
use MOM_error_handler, only : MOM_error, FATAL, WARNING, callTree_showQuery
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,  only : forcing, extractFluxes1d, forcing_SinglePointPrint
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : slasher
use MOM_opacity,       only : set_opacity, opacity_CS, extract_optics_slice, extract_optics_fields
use MOM_opacity,       only : optics_type, optics_nbands, absorbRemainingSW, sumSWoverBands
use MOM_tracer_flow_control, only : get_chl_from_model, tracer_flow_control_CS
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs ! , vertvisc_type, accel_diag_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use time_interp_external_mod, only : init_external_field, time_interp_external
use time_interp_external_mod, only : time_interp_external_init

implicit none ; private

#include <MOM_memory.h>

public diabatic_aux_init, diabatic_aux_end
public make_frazil, adjust_salt, differential_diffuse_T_S, triDiagTS, triDiagTS_Eulerian
public find_uv_at_h, diagnoseMLDbyDensityDifference, applyBoundaryFluxesInOut, set_pen_shortwave
public diagnoseMLDbyEnergy

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for diabatic_aux
type, public :: diabatic_aux_CS ; private
  logical :: do_rivermix = .false. !< Provide additional TKE to mix river runoff at the
                                   !! river mouths to a depth of "rivermix_depth"
  real    :: rivermix_depth = 0.0  !< The depth to which rivers are mixed if do_rivermix = T [Z ~> m].
  logical :: reclaim_frazil  !<   If true, try to use any frazil heat deficit to
                             !! to cool the topmost layer down to the freezing
                             !! point.  The default is false.
  logical :: pressure_dependent_frazil  !< If true, use a pressure dependent
                             !! freezing temperature when making frazil.  The
                             !! default is false, which will be faster but is
                             !! inappropriate with ice-shelf cavities.
  logical :: ignore_fluxes_over_land    !< If true, the model does not check
                             !! if fluxes are applied over land points. This
                             !! flag must be used when the ocean is coupled with
                             !! sea ice and ice shelves and use_ePBL = true.
  logical :: use_river_heat_content !< If true, assumes that ice-ocean boundary
                             !! has provided a river heat content. Otherwise, runoff
                             !! is added with a temperature of the local SST.
  logical :: use_calving_heat_content !< If true, assumes that ice-ocean boundary
                             !! has provided a calving heat content. Otherwise, calving
                             !! is added with a temperature of the local SST.
  logical :: var_pen_sw      !<   If true, use one of the CHL_A schemes to determine the
                             !! e-folding depth of incoming shortwave radiation.
  integer :: sbc_chl         !< An integer handle used in time interpolation of
                             !! chlorophyll read from a file.
  logical ::  chl_from_file  !< If true, chl_a is read from a file.

  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< Structure used to regulate timing of diagnostic output

  ! Diagnostic handles
  integer :: id_createdH       = -1 !< Diagnostic ID of mass added to avoid grounding
  integer :: id_brine_lay      = -1 !< Diagnostic ID of which layer receives the brine
  integer :: id_penSW_diag     = -1 !< Diagnostic ID of Penetrative shortwave heating (flux convergence)
  integer :: id_penSWflux_diag = -1 !< Diagnostic ID of Penetrative shortwave flux
  integer :: id_nonpenSW_diag  = -1 !< Diagnostic ID of Non-penetrative shortwave heating
  integer :: id_Chl            = -1 !< Diagnostic ID of chlorophyll-A handles for opacity

  ! Optional diagnostic arrays
  real, allocatable, dimension(:,:)   :: createdH       !< The amount of volume added in order to
                                                        !! avoid grounding [H T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:) :: penSW_diag     !< Heating in a layer from convergence of
                                                        !! penetrative SW [Q R Z T-1 ~> W m-2]
  real, allocatable, dimension(:,:,:) :: penSWflux_diag !< Penetrative SW flux at base of grid
                                                        !! layer [Q R Z T-1 ~> W m-2]
  real, allocatable, dimension(:,:)   :: nonpenSW_diag  !< Non-downwelling SW radiation at ocean
                                                        !! surface [Q R Z T-1 ~> W m-2]

end type diabatic_aux_CS

!>@{ CPU time clock IDs
integer :: id_clock_uv_at_h, id_clock_frazil
!>@}

contains

!> Frazil formation keeps the temperature above the freezing point.
!! This subroutine warms any water that is colder than the (currently
!! surface) freezing point up to the freezing point and accumulates
!! the required heat (in [Q R Z ~> J m-2]) in tv%frazil.
subroutine make_frazil(h, tv, G, GV, US, CS, p_surf, halo)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv !< Structure containing pointers to any available
                                               !! thermodynamic fields.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(diabatic_aux_CS),   intent(in)    :: CS !< The control structure returned by a previous
                                               !! call to diabatic_aux_init.
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa].
  integer,       optional, intent(in)    :: halo !< Halo width over which to calculate frazil

  ! Local variables
  real, dimension(SZI_(G)) :: &
    fraz_col, & ! The accumulated heat requirement due to frazil [Q R Z ~> J m-2].
    T_freeze, & ! The freezing potential temperature at the current salinity [degC].
    ps          ! Surface pressure [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G),SZK_(G)) :: &
    pressure    ! The pressure at the middle of each layer [R L2 T-2 ~> Pa].
  real :: H_to_RL2_T2  ! A conversion factor from thicknesses in H to pressure [R L2 T-2 H-1 ~> Pa m-1 or Pa m2 kg-1]
  real :: hc    ! A layer's heat capacity [Q R Z degC-1 ~> J m-2 degC-1].
  logical :: T_fr_set  ! True if the freezing point has been calculated for a
                       ! row of points.
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  endif

  call cpu_clock_begin(id_clock_frazil)

  if (.not.CS%pressure_dependent_frazil) then
    do k=1,nz ; do i=is,ie ; pressure(i,k) = 0.0 ; enddo ; enddo
  else
    H_to_RL2_T2 = GV%H_to_RZ * GV%g_Earth
  endif
  !$OMP parallel do default(shared) private(fraz_col,T_fr_set,T_freeze,hc,ps)  &
  !$OMP                             firstprivate(pressure) ! pressure might be set above, so should be firstprivate
  do j=js,je
    ps(:) = 0.0
    if (PRESENT(p_surf)) then ; do i=is,ie
      ps(i) = p_surf(i,j)
    enddo ; endif

    do i=is,ie ; fraz_col(i) = 0.0 ; enddo

    if (CS%pressure_dependent_frazil) then
      do i=is,ie
        pressure(i,1) = ps(i) + (0.5*H_to_RL2_T2)*h(i,j,1)
      enddo
      do k=2,nz ; do i=is,ie
        pressure(i,k) = pressure(i,k-1) + &
          (0.5*H_to_RL2_T2) * (h(i,j,k) + h(i,j,k-1))
      enddo ; enddo
    endif

    if (CS%reclaim_frazil) then
      T_fr_set = .false.
      do i=is,ie ; if (tv%frazil(i,j) > 0.0) then
        if (.not.T_fr_set) then
          call calculate_TFreeze(tv%S(i:,j,1), pressure(i:,1), T_freeze(i:), &
                                 1, ie-i+1, tv%eqn_of_state, pres_scale=US%RL2_T2_to_Pa)
          T_fr_set = .true.
        endif

        if (tv%T(i,j,1) > T_freeze(i)) then
    ! If frazil had previously been formed, but the surface temperature is now
    ! above freezing, cool the surface layer with the frazil heat deficit.
          hc = (tv%C_p*GV%H_to_RZ) * h(i,j,1)
          if (tv%frazil(i,j) - hc * (tv%T(i,j,1) - T_freeze(i)) <= 0.0) then
            tv%T(i,j,1) = tv%T(i,j,1) - tv%frazil(i,j) / hc
            tv%frazil(i,j) = 0.0
          else
            tv%frazil(i,j) = tv%frazil(i,j) - hc * (tv%T(i,j,1) - T_freeze(i))
            tv%T(i,j,1) = T_freeze(i)
          endif
        endif
      endif ; enddo
    endif

    do k=nz,1,-1
      T_fr_set = .false.
      do i=is,ie
        if ((G%mask2dT(i,j) > 0.0) .and. &
            ((tv%T(i,j,k) < 0.0) .or. (fraz_col(i) > 0.0))) then
          if (.not.T_fr_set) then
            call calculate_TFreeze(tv%S(i:,j,k), pressure(i:,k), T_freeze(i:), &
                                   1, ie-i+1, tv%eqn_of_state, pres_scale=US%RL2_T2_to_Pa)
            T_fr_set = .true.
          endif

          hc = (tv%C_p*GV%H_to_RZ) * h(i,j,k)
          if (h(i,j,k) <= 10.0*(GV%Angstrom_H + GV%H_subroundoff)) then
            ! Very thin layers should not be cooled by the frazil flux.
            if (tv%T(i,j,k) < T_freeze(i)) then
              fraz_col(i) = fraz_col(i) + hc * (T_freeze(i) - tv%T(i,j,k))
              tv%T(i,j,k) = T_freeze(i)
            endif
          elseif ((fraz_col(i) > 0.0) .or. (tv%T(i,j,k) < T_freeze(i))) then
            if (fraz_col(i) + hc * (T_freeze(i) - tv%T(i,j,k)) < 0.0) then
              tv%T(i,j,k) = tv%T(i,j,k) - fraz_col(i) / hc
              fraz_col(i) = 0.0
            else
              fraz_col(i) = fraz_col(i) + hc * (T_freeze(i) - tv%T(i,j,k))
              tv%T(i,j,k) = T_freeze(i)
            endif
          endif
        endif
      enddo
    enddo
    do i=is,ie
      tv%frazil(i,j) = tv%frazil(i,j) + fraz_col(i)
    enddo
  enddo
  call cpu_clock_end(id_clock_frazil)

end subroutine make_frazil

!> This subroutine applies double diffusion to T & S, assuming no diapycal mass
!! fluxes, using a simple triadiagonal solver.
subroutine differential_diffuse_T_S(h, T, S, Kd_T, Kd_S, dt, G, GV)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: T    !< Potential temperature [degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: S    !< Salinity [PSU] or [gSalt/kg], generically [ppt].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), &
                           intent(inout)    :: Kd_T !< The extra diffusivity of temperature due to
                                                 !! double diffusion relative to the diffusivity of
                                                 !! diffusivity of density [Z2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), &
                           intent(in)    :: Kd_S !< The extra diffusivity of salinity due to
                                                 !! double diffusion relative to the diffusivity of
                                                 !! diffusivity of density [Z2 T-1 ~> m2 s-1].
  real,                    intent(in)    :: dt   !<  Time increment [T ~> s].

  ! local variables
  real, dimension(SZI_(G)) :: &
    b1_T, b1_S, &  !  Variables used by the tridiagonal solvers of T & S [H ~> m or kg m-2].
    d1_T, d1_S     !  Variables used by the tridiagonal solvers [nondim].
  real, dimension(SZI_(G),SZK_(G)) :: &
    c1_T, c1_S     !  Variables used by the tridiagonal solvers [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(G)+1) :: &
    mix_T, mix_S   !  Mixing distances in both directions across each interface [H ~> m or kg m-2].
  real :: h_tr         ! h_tr is h at tracer points with a tiny thickness
                       ! added to ensure positive definiteness [H ~> m or kg m-2].
  real :: h_neglect    ! A thickness that is so small it is usually lost
                       ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: I_h_int      ! The inverse of the thickness associated with an
                       ! interface [H-1 ~> m-1 or m2 kg-1].
  real :: b_denom_T    ! The first term in the denominators for the expressions
  real :: b_denom_S    ! for b1_T and b1_S, both [H ~> m or kg m-2].
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  h_neglect = GV%H_subroundoff

  !$OMP parallel do default(private) shared(is,ie,js,je,h,h_neglect,dt,Kd_T,Kd_S,G,GV,T,S,nz)
  do j=js,je
    do i=is,ie
      I_h_int = 1.0 / (0.5 * (h(i,j,1) + h(i,j,2)) + h_neglect)
      mix_T(i,2) = ((dt * Kd_T(i,j,2)) * GV%Z_to_H**2) * I_h_int
      mix_S(i,2) = ((dt * Kd_S(i,j,2)) * GV%Z_to_H**2) * I_h_int

      h_tr = h(i,j,1) + h_neglect
      b1_T(i) = 1.0 / (h_tr + mix_T(i,2))
      b1_S(i) = 1.0 / (h_tr + mix_S(i,2))
      d1_T(i) = h_tr * b1_T(i)
      d1_S(i) = h_tr * b1_S(i)
      T(i,j,1) = (b1_T(i)*h_tr)*T(i,j,1)
      S(i,j,1) = (b1_S(i)*h_tr)*S(i,j,1)
    enddo
    do k=2,nz-1 ; do i=is,ie
      ! Calculate the mixing across the interface below this layer.
      I_h_int = 1.0 / (0.5 * (h(i,j,k) + h(i,j,k+1)) + h_neglect)
      mix_T(i,K+1) = ((dt * Kd_T(i,j,K+1)) * GV%Z_to_H**2) * I_h_int
      mix_S(i,K+1) = ((dt * Kd_S(i,j,K+1)) * GV%Z_to_H**2) * I_h_int

      c1_T(i,k) = mix_T(i,K) * b1_T(i)
      c1_S(i,k) = mix_S(i,K) * b1_S(i)

      h_tr = h(i,j,k) + h_neglect
      b_denom_T = h_tr + d1_T(i)*mix_T(i,K)
      b_denom_S = h_tr + d1_S(i)*mix_S(i,K)
      b1_T(i) = 1.0 / (b_denom_T + mix_T(i,K+1))
      b1_S(i) = 1.0 / (b_denom_S + mix_S(i,K+1))
      d1_T(i) = b_denom_T * b1_T(i)
      d1_S(i) = b_denom_S * b1_S(i)

      T(i,j,k) = b1_T(i) * (h_tr*T(i,j,k) + mix_T(i,K)*T(i,j,k-1))
      S(i,j,k) = b1_S(i) * (h_tr*S(i,j,k) + mix_S(i,K)*S(i,j,k-1))
    enddo ; enddo
    do i=is,ie
      c1_T(i,nz) = mix_T(i,nz) * b1_T(i)
      c1_S(i,nz) = mix_S(i,nz) * b1_S(i)

      h_tr = h(i,j,nz) + h_neglect
      b1_T(i) = 1.0 / (h_tr + d1_T(i)*mix_T(i,nz))
      b1_S(i) = 1.0 / (h_tr + d1_S(i)*mix_S(i,nz))

      T(i,j,nz) = b1_T(i) * (h_tr*T(i,j,nz) + mix_T(i,nz)*T(i,j,nz-1))
      S(i,j,nz) = b1_S(i) * (h_tr*S(i,j,nz) + mix_S(i,nz)*S(i,j,nz-1))
    enddo
    do k=nz-1,1,-1 ; do i=is,ie
      T(i,j,k) = T(i,j,k) + c1_T(i,k+1)*T(i,j,k+1)
      S(i,j,k) = S(i,j,k) + c1_S(i,k+1)*S(i,j,k+1)
    enddo ; enddo
  enddo
end subroutine differential_diffuse_T_S

!> This subroutine keeps salinity from falling below a small but positive threshold.
!! This usually occurs when the ice model attempts to extract more salt then
!! is actually available to it from the ocean.
subroutine adjust_salt(h, tv, G, GV, CS, halo)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv   !< Structure containing pointers to any
                                                 !! available thermodynamic fields.
  type(diabatic_aux_CS),   intent(in)    :: CS   !< The control structure returned by a previous
                                                 !! call to diabatic_aux_init.
  integer,       optional, intent(in)    :: halo !< Halo width over which to work

  ! local variables
  real :: salt_add_col(SZI_(G),SZJ_(G)) !< The accumulated salt requirement [ppt R Z ~> gSalt m-2]
  real :: S_min      !< The minimum salinity [ppt].
  real :: mc         !< A layer's mass [R Z ~> kg m-2].
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  endif

!  call cpu_clock_begin(id_clock_adjust_salt)

  S_min = tv%min_salinity

  salt_add_col(:,:) = 0.0

  !$OMP parallel do default(shared) private(mc)
  do j=js,je
    do k=nz,1,-1 ; do i=is,ie
      if ( (G%mask2dT(i,j) > 0.0) .and. &
           ((tv%S(i,j,k) < S_min) .or. (salt_add_col(i,j) > 0.0)) ) then
        mc = GV%H_to_RZ * h(i,j,k)
        if (h(i,j,k) <= 10.0*GV%Angstrom_H) then
          ! Very thin layers should not be adjusted by the salt flux
          if (tv%S(i,j,k) < S_min) then
            salt_add_col(i,j) = salt_add_col(i,j) +  mc * (S_min - tv%S(i,j,k))
            tv%S(i,j,k) = S_min
          endif
        elseif (salt_add_col(i,j) + mc * (S_min - tv%S(i,j,k)) <= 0.0) then
          tv%S(i,j,k) = tv%S(i,j,k) - salt_add_col(i,j) / mc
          salt_add_col(i,j) = 0.0
        else
          salt_add_col(i,j) = salt_add_col(i,j) + mc * (S_min - tv%S(i,j,k))
          tv%S(i,j,k) = S_min
        endif
      endif
    enddo ; enddo
    do i=is,ie
      tv%salt_deficit(i,j) = tv%salt_deficit(i,j) + salt_add_col(i,j)
    enddo
  enddo
!  call cpu_clock_end(id_clock_adjust_salt)

end subroutine adjust_salt

!> This is a simple tri-diagonal solver for T and S.
!! "Simple" means it only uses arrays hold, ea and eb.
subroutine triDiagTS(G, GV, is, ie, js, je, hold, ea, eb, T, S)
  type(ocean_grid_type),                    intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                  intent(in)    :: GV   !< The ocean's vertical grid structure
  integer,                                  intent(in)    :: is   !< The start i-index to work on.
  integer,                                  intent(in)    :: ie   !< The end i-index to work on.
  integer,                                  intent(in)    :: js   !< The start j-index to work on.
  integer,                                  intent(in)    :: je   !< The end j-index to work on.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: hold !< The layer thicknesses before entrainment,
                                                                  !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: ea !< The amount of fluid entrained from the layer
                                                                !! above within this time step [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: eb !< The amount of fluid entrained from the layer
                                                                !! below within this time step [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: T  !< Layer potential temperatures [degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: S  !< Layer salinities [ppt].

  ! Local variables
  real :: b1(SZIB_(G))          ! A variable used by the tridiagonal solver [H-1 ~> m-2 or m2 kg-1].
  real :: d1(SZIB_(G))          ! A variable used by the tridiagonal solver [nondim].
  real :: c1(SZIB_(G),SZK_(G))  ! A variable used by the tridiagonal solver [nondim].
  real :: h_tr, b_denom_1       ! Two temporary thicknesses [H ~> m or kg m-2].
  integer :: i, j, k

  !$OMP parallel do default(shared) private(h_tr,b1,d1,c1,b_denom_1)
  do j=js,je
    do i=is,ie
      h_tr = hold(i,j,1) + GV%H_subroundoff
      b1(i) = 1.0 / (h_tr + eb(i,j,1))
      d1(i) = h_tr * b1(i)
      T(i,j,1) = (b1(i)*h_tr)*T(i,j,1)
      S(i,j,1) = (b1(i)*h_tr)*S(i,j,1)
    enddo
    do k=2,G%ke ; do i=is,ie
      c1(i,k) = eb(i,j,k-1) * b1(i)
      h_tr = hold(i,j,k) + GV%H_subroundoff
      b_denom_1 = h_tr + d1(i)*ea(i,j,k)
      b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
      d1(i) = b_denom_1 * b1(i)
      T(i,j,k) = b1(i) * (h_tr*T(i,j,k) + ea(i,j,k)*T(i,j,k-1))
      S(i,j,k) = b1(i) * (h_tr*S(i,j,k) + ea(i,j,k)*S(i,j,k-1))
    enddo ; enddo
    do k=G%ke-1,1,-1 ; do i=is,ie
      T(i,j,k) = T(i,j,k) + c1(i,k+1)*T(i,j,k+1)
      S(i,j,k) = S(i,j,k) + c1(i,k+1)*S(i,j,k+1)
    enddo ; enddo
  enddo
end subroutine triDiagTS

!> This is a simple tri-diagonal solver for T and S, with mixing across interfaces but no net
!! transfer of mass.
subroutine triDiagTS_Eulerian(G, GV, is, ie, js, je, hold, ent, T, S)
  type(ocean_grid_type),                    intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                  intent(in)    :: GV   !< The ocean's vertical grid structure
  integer,                                  intent(in)    :: is   !< The start i-index to work on.
  integer,                                  intent(in)    :: ie   !< The end i-index to work on.
  integer,                                  intent(in)    :: js   !< The start j-index to work on.
  integer,                                  intent(in)    :: je   !< The end j-index to work on.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: hold !< The layer thicknesses before entrainment,
                                                                  !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(in)  :: ent  !< The amount of fluid mixed across an interface
                                                                  !! within this time step [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: T    !< Layer potential temperatures [degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: S    !< Layer salinities [ppt].

  ! Local variables
  real :: b1(SZIB_(G))          ! A variable used by the tridiagonal solver [H-1 ~> m-2 or m2 kg-1].
  real :: d1(SZIB_(G))          ! A variable used by the tridiagonal solver [nondim].
  real :: c1(SZIB_(G),SZK_(G))  ! A variable used by the tridiagonal solver [nondim].
  real :: h_tr, b_denom_1       ! Two temporary thicknesses [H ~> m or kg m-2].
  integer :: i, j, k

  !$OMP parallel do default(shared) private(h_tr,b1,d1,c1,b_denom_1)
  do j=js,je
    do i=is,ie
      h_tr = hold(i,j,1) + GV%H_subroundoff
      b1(i) = 1.0 / (h_tr + ent(i,j,2))
      d1(i) = h_tr * b1(i)
      T(i,j,1) = (b1(i)*h_tr)*T(i,j,1)
      S(i,j,1) = (b1(i)*h_tr)*S(i,j,1)
    enddo
    do k=2,G%ke ; do i=is,ie
      c1(i,k) = ent(i,j,K) * b1(i)
      h_tr = hold(i,j,k) + GV%H_subroundoff
      b_denom_1 = h_tr + d1(i)*ent(i,j,K)
      b1(i) = 1.0 / (b_denom_1 + ent(i,j,K+1))
      d1(i) = b_denom_1 * b1(i)
      T(i,j,k) = b1(i) * (h_tr*T(i,j,k) + ent(i,j,K)*T(i,j,k-1))
      S(i,j,k) = b1(i) * (h_tr*S(i,j,k) + ent(i,j,K)*S(i,j,k-1))
    enddo ; enddo
    do k=G%ke-1,1,-1 ; do i=is,ie
      T(i,j,k) = T(i,j,k) + c1(i,k+1)*T(i,j,k+1)
      S(i,j,k) = S(i,j,k) + c1(i,k+1)*S(i,j,k+1)
    enddo ; enddo
  enddo
end subroutine triDiagTS_Eulerian


!>   This subroutine calculates u_h and v_h (velocities at thickness
!! points), optionally using the entrainment amounts passed in as arguments.
subroutine find_uv_at_h(u, v, h, u_h, v_h, G, GV, US, ea, eb, zero_mix)
  type(ocean_grid_type),     intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),   intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),     intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                             intent(in)  :: u    !< The zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                             intent(in)  :: v    !< The meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                             intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                             intent(out)   :: u_h !< Zonal velocity interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                             intent(out)   :: v_h !< Meridional velocity interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                     optional, intent(in)  :: ea !< The amount of fluid entrained from the layer
                                                 !! above within this time step [H ~> m or kg m-2].
                                                 !! Omitting ea is the same as setting it to 0.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                     optional, intent(in)  :: eb !< The amount of fluid entrained from the layer
                                                 !! below within this time step [H ~> m or kg m-2].
                                                 !! Omitting eb is the same as setting it to 0.
  logical,           optional, intent(in)  :: zero_mix !< If true, do the calculation of u_h and
                                                 !! v_h as though ea and eb were being supplied with
                                                 !! uniformly zero values.

  ! local variables
  real :: b_denom_1    ! The first term in the denominator of b1 [H ~> m or kg m-2].
  real :: h_neglect    ! A thickness that is so small it is usually lost
                       ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: b1(SZI_(G)), d1(SZI_(G)), c1(SZI_(G),SZK_(G))
  real :: a_n(SZI_(G)), a_s(SZI_(G))  ! Fractional weights of the neighboring
  real :: a_e(SZI_(G)), a_w(SZI_(G))  ! velocity points, ~1/2 in the open
                                      ! ocean, nondimensional.
  real :: sum_area, Idenom
  logical :: mix_vertically, zero_mixing
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  call cpu_clock_begin(id_clock_uv_at_h)
  h_neglect = GV%H_subroundoff

  mix_vertically = present(ea)
  if (present(ea) .neqv. present(eb)) call MOM_error(FATAL, &
      "find_uv_at_h: Either both ea and eb or neither one must be present "// &
      "in call to find_uv_at_h.")
  zero_mixing = .false. ; if (present(zero_mix)) zero_mixing = zero_mix
  if (zero_mixing) mix_vertically = .false.
  !$OMP parallel do default(none) shared(is,ie,js,je,G,GV,mix_vertically,zero_mixing,h, &
  !$OMP                                  h_neglect,ea,eb,u_h,u,v_h,v,nz)                &
  !$OMP                          private(sum_area,Idenom,a_w,a_e,a_s,a_n,b_denom_1,b1,d1,c1)
  do j=js,je
    do i=is,ie
      sum_area = G%areaCu(I-1,j) + G%areaCu(I,j)
      if (sum_area>0.0) then
        Idenom = sqrt(0.5*G%IareaT(i,j) / sum_area)
        a_w(i) = G%areaCu(I-1,j) * Idenom
        a_e(i) = G%areaCu(I,j) * Idenom
      else
        a_w(i) = 0.0 ; a_e(i) = 0.0
      endif

      sum_area = G%areaCv(i,J-1) + G%areaCv(i,J)
      if (sum_area>0.0) then
        Idenom = sqrt(0.5*G%IareaT(i,j) / sum_area)
        a_s(i) = G%areaCv(i,J-1) * Idenom
        a_n(i) = G%areaCv(i,J) * Idenom
      else
        a_s(i) = 0.0 ; a_n(i) = 0.0
      endif
    enddo

    if (mix_vertically) then
      do i=is,ie
        b_denom_1 = h(i,j,1) + h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,1))
        d1(i) = b_denom_1 * b1(i)
        u_h(i,j,1) = (h(i,j,1)*b1(i)) * (a_e(i)*u(I,j,1) + a_w(i)*u(I-1,j,1))
        v_h(i,j,1) = (h(i,j,1)*b1(i)) * (a_n(i)*v(i,J,1) + a_s(i)*v(i,J-1,1))
      enddo
      do k=2,nz ; do i=is,ie
        c1(i,k) = eb(i,j,k-1) * b1(i)
        b_denom_1 = h(i,j,k) + d1(i)*ea(i,j,k) + h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
        d1(i) = b_denom_1 * b1(i)
        u_h(i,j,k) = (h(i,j,k) * (a_e(i)*u(I,j,k) + a_w(i)*u(I-1,j,k)) + &
                      ea(i,j,k)*u_h(i,j,k-1))*b1(i)
        v_h(i,j,k) = (h(i,j,k) * (a_n(i)*v(i,J,k) + a_s(i)*v(i,J-1,k)) + &
                      ea(i,j,k)*v_h(i,j,k-1))*b1(i)
      enddo ; enddo
      do k=nz-1,1,-1 ; do i=is,ie
        u_h(i,j,k) = u_h(i,j,k) + c1(i,k+1)*u_h(i,j,k+1)
        v_h(i,j,k) = v_h(i,j,k) + c1(i,k+1)*v_h(i,j,k+1)
      enddo ; enddo
    elseif (zero_mixing) then
      do i=is,ie
        b1(i) = 1.0 / (h(i,j,1) + h_neglect)
        u_h(i,j,1) = (h(i,j,1)*b1(i)) * (a_e(i)*u(I,j,1) + a_w(i)*u(I-1,j,1))
        v_h(i,j,1) = (h(i,j,1)*b1(i)) * (a_n(i)*v(i,J,1) + a_s(i)*v(i,J-1,1))
      enddo
      do k=2,nz ; do i=is,ie
        b1(i) = 1.0 / (h(i,j,k) + h_neglect)
        u_h(i,j,k) = (h(i,j,k) * (a_e(i)*u(I,j,k) + a_w(i)*u(I-1,j,k))) * b1(i)
        v_h(i,j,k) = (h(i,j,k) * (a_n(i)*v(i,J,k) + a_s(i)*v(i,J-1,k))) * b1(i)
      enddo ; enddo
    else
      do k=1,nz ; do i=is,ie
        u_h(i,j,k) = a_e(i)*u(I,j,k) + a_w(i)*u(I-1,j,k)
        v_h(i,j,k) = a_n(i)*v(i,J,k) + a_s(i)*v(i,J-1,k)
      enddo ; enddo
    endif
  enddo

  call cpu_clock_end(id_clock_uv_at_h)
end subroutine find_uv_at_h


subroutine set_pen_shortwave(optics, fluxes, G, GV, US, CS, opacity_CSp, tracer_flow_CSp)
  type(optics_type),       pointer       :: optics !< An optics structure that has will contain
                                                   !! information about shortwave fluxes and absorption.
  type(forcing),           intent(inout) :: fluxes !< points to forcing fields
                                                   !! unused fields have NULL ptrs
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(diabatic_aux_CS),   pointer       :: CS     !< Control structure for diabatic_aux
  type(opacity_CS),        pointer       :: opacity_CSp !< The control structure for the opacity module.
  type(tracer_flow_control_CS), pointer  :: tracer_flow_CSp !< A pointer to the control structure
                                                   !! organizing the tracer modules.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))          :: chl_2d !< Vertically uniform chlorophyll-A concentractions [mg m-3]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: chl_3d !< The chlorophyll-A concentractions of each layer [mg m-3]
  character(len=128) :: mesg
  integer :: i, j, k, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.not.associated(optics)) return

  if (CS%var_pen_sw) then
    if (CS%chl_from_file) then
      ! Only the 2-d surface chlorophyll can be read in from a file.  The
      ! same value is assumed for all layers.
      call time_interp_external(CS%sbc_chl, CS%Time, chl_2d)
      do j=js,je ; do i=is,ie
        if ((G%mask2dT(i,j) > 0.5) .and. (chl_2d(i,j) < 0.0)) then
          write(mesg,'(" Time_interp negative chl of ",(1pe12.4)," at i,j = ",&
                    & 2(i3), "lon/lat = ",(1pe12.4)," E ", (1pe12.4), " N.")') &
                     chl_2d(i,j), i, j, G%geoLonT(i,j), G%geoLatT(i,j)
          call MOM_error(FATAL, "MOM_diabatic_aux set_pen_shortwave: "//trim(mesg))
        endif
      enddo ; enddo

      if (CS%id_chl > 0) call post_data(CS%id_chl, chl_2d, CS%diag)

      call set_opacity(optics, fluxes%sw, fluxes%sw_vis_dir, fluxes%sw_vis_dif, &
                       fluxes%sw_nir_dir, fluxes%sw_nir_dif, G, GV, US, opacity_CSp, chl_2d=chl_2d)
    else
      if (.not.associated(tracer_flow_CSp)) call MOM_error(FATAL, &
        "The tracer flow control structure must be associated when the model sets "//&
        "the chlorophyll internally in set_pen_shortwave.")
      call get_chl_from_model(chl_3d, G, tracer_flow_CSp)

      if (CS%id_chl > 0) call post_data(CS%id_chl, chl_3d(:,:,1), CS%diag)

      call set_opacity(optics, fluxes%sw, fluxes%sw_vis_dir, fluxes%sw_vis_dif, &
                       fluxes%sw_nir_dir, fluxes%sw_nir_dif, G, GV, US, opacity_CSp, chl_3d=chl_3d)
    endif
  else
    call set_opacity(optics, fluxes%sw, fluxes%sw_vis_dir, fluxes%sw_vis_dif, &
                     fluxes%sw_nir_dir, fluxes%sw_nir_dif, G, GV, US, opacity_CSp)
  endif

end subroutine set_pen_shortwave


!> Diagnose a mixed layer depth (MLD) determined by a given density difference with the surface.
!> This routine is appropriate in MOM_diabatic_driver due to its position within the time stepping.
subroutine diagnoseMLDbyDensityDifference(id_MLD, h, tv, densityDiff, G, GV, US, diagPtr, &
                                          id_N2subML, id_MLDsq, dz_subML)
  type(ocean_grid_type),   intent(in) :: G           !< Grid type
  type(verticalGrid_type), intent(in) :: GV          !< ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US          !< A dimensional unit scaling type
  integer,                 intent(in) :: id_MLD      !< Handle (ID) of MLD diagnostic
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h           !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in) :: tv          !< Structure containing pointers to any
                                                     !! available thermodynamic fields.
  real,                    intent(in) :: densityDiff !< Density difference to determine MLD [R ~> kg m-3]
  type(diag_ctrl),         pointer    :: diagPtr     !< Diagnostics structure
  integer,       optional, intent(in) :: id_N2subML  !< Optional handle (ID) of subML stratification
  integer,       optional, intent(in) :: id_MLDsq    !< Optional handle (ID) of squared MLD
  real,          optional, intent(in) :: dz_subML    !< The distance over which to calculate N2subML
                                                     !! or 50 m if missing [Z ~> m]

  ! Local variables
  real, dimension(SZI_(G)) :: deltaRhoAtKm1, deltaRhoAtK ! Density differences [R ~> kg m-3].
  real, dimension(SZI_(G)) :: pRef_MLD, pRef_N2 ! Reference pressures [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G)) :: H_subML, dH_N2    ! Summed thicknesses used in N2 calculation [H ~> m].
  real, dimension(SZI_(G)) :: T_subML, T_deeper ! Temperatures used in the N2 calculation [degC].
  real, dimension(SZI_(G)) :: S_subML, S_deeper ! Salinities used in the N2 calculation [ppt].
  real, dimension(SZI_(G)) :: rho_subML, rho_deeper ! Densities used in the N2 calculation [R ~> kg m-3].
  real, dimension(SZI_(G)) :: dK, dKm1          ! Depths [Z ~> m].
  real, dimension(SZI_(G)) :: rhoSurf          ! Density used in finding the mixedlayer depth [R ~> kg m-3].
  real, dimension(SZI_(G), SZJ_(G)) :: MLD     ! Diagnosed mixed layer depth [Z ~> m].
  real, dimension(SZI_(G), SZJ_(G)) :: subMLN2 ! Diagnosed stratification below ML [T-2 ~> s-2].
  real, dimension(SZI_(G), SZJ_(G)) :: MLD2    ! Diagnosed MLD^2 [Z2 ~> m2].
  logical, dimension(SZI_(G)) :: N2_region_set ! If true, all necessary values for calculating N2
                                               ! have been stored already.
  real :: gE_Rho0          ! The gravitational acceleration divided by a mean density [Z T-2 R-1 ~> m4 s-2 kg-1].
  real :: dH_subML         ! Depth below ML over which to diagnose stratification [H ~> m].
  real :: aFac             ! A nondimensional factor [nondim]
  real :: ddRho            ! A density difference [R ~> kg m-3]
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, is, ie, js, je, k, nz, id_N2, id_SQ

  id_N2 = -1 ; if (PRESENT(id_N2subML)) id_N2 = id_N2subML

  id_SQ = -1 ; if (PRESENT(id_MLDsq)) id_SQ = id_MLDsq

  gE_rho0 = US%L_to_Z**2*GV%g_Earth / (GV%Rho0)
  dH_subML = 50.*GV%m_to_H  ; if (present(dz_subML)) dH_subML = GV%Z_to_H*dz_subML

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  pRef_MLD(:) = 0.0
  EOSdom(:) = EOS_domain(G%HI)
  do j=js,je
    do i=is,ie ; dK(i) = 0.5 * h(i,j,1) * GV%H_to_Z ; enddo ! Depth of center of surface layer
    call calculate_density(tv%T(:,j,1), tv%S(:,j,1), pRef_MLD, rhoSurf, tv%eqn_of_state, EOSdom)
    do i=is,ie
      deltaRhoAtK(i) = 0.
      MLD(i,j) = 0.
      if (id_N2>0) then
        subMLN2(i,j) = 0.0
        H_subML(i) = h(i,j,1) ; dH_N2(i) = 0.0
        T_subML(i) = 0.0  ; S_subML(i) = 0.0 ; T_deeper(i) = 0.0 ; S_deeper(i) = 0.0
        N2_region_set(i) = (G%mask2dT(i,j)<0.5) ! Only need to work on ocean points.
      endif
    enddo
    do k=2,nz
      do i=is,ie
        dKm1(i) = dK(i) ! Depth of center of layer K-1
        dK(i) = dK(i) + 0.5 * ( h(i,j,k) + h(i,j,k-1) ) * GV%H_to_Z ! Depth of center of layer K
      enddo

      ! Prepare to calculate stratification, N2, immediately below the mixed layer by finding
      ! the cells that extend over at least dz_subML.
      if (id_N2>0) then
         do i=is,ie
          if (MLD(i,j)==0.0) then  ! Still in the mixed layer.
            H_subML(i) = H_subML(i) + h(i,j,k)
          elseif (.not.N2_region_set(i)) then ! This block is below the mixed layer, but N2 has not been found yet.
            if (dH_N2(i)==0.0) then ! Record the temperature, salinity, pressure, immediately below the ML
              T_subML(i) = tv%T(i,j,k) ; S_subML(i) = tv%S(i,j,k)
              H_subML(i) = H_subML(i) + 0.5 * h(i,j,k) ! Start midway through this layer.
              dH_N2(i) = 0.5 * h(i,j,k)
            elseif (dH_N2(i) + h(i,j,k) < dH_subML) then
              dH_N2(i) = dH_N2(i) + h(i,j,k)
            else  ! This layer includes the base of the region where N2 is calculated.
              T_deeper(i) = tv%T(i,j,k) ; S_deeper(i) = tv%S(i,j,k)
              dH_N2(i) = dH_N2(i) + 0.5 * h(i,j,k)
              N2_region_set(i) = .true.
            endif
          endif
        enddo ! i-loop
      endif ! id_N2>0

      ! Mixed-layer depth, using sigma-0 (surface reference pressure)
      do i=is,ie ; deltaRhoAtKm1(i) = deltaRhoAtK(i) ; enddo ! Store value from previous iteration of K
      call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pRef_MLD, deltaRhoAtK, tv%eqn_of_state, EOSdom)
      do i = is, ie
        deltaRhoAtK(i) = deltaRhoAtK(i) - rhoSurf(i) ! Density difference between layer K and surface
        ddRho = deltaRhoAtK(i) - deltaRhoAtKm1(i)
        if ((MLD(i,j)==0.) .and. (ddRho>0.) .and. &
            (deltaRhoAtKm1(i)<densityDiff) .and. (deltaRhoAtK(i)>=densityDiff)) then
          aFac = ( densityDiff - deltaRhoAtKm1(i) ) / ddRho
          MLD(i,j) = dK(i) * aFac + dKm1(i) * (1. - aFac)
        endif
        if (id_SQ > 0) MLD2(i,j) = MLD(i,j)**2
      enddo ! i-loop
    enddo ! k-loop
    do i=is,ie
      if ((MLD(i,j)==0.) .and. (deltaRhoAtK(i)<densityDiff)) MLD(i,j) = dK(i) ! Assume mixing to the bottom
    enddo

    if (id_N2>0) then  ! Now actually calculate stratification, N2, below the mixed layer.
      do i=is,ie ; pRef_N2(i) = (GV%g_Earth * GV%H_to_RZ) * (H_subML(i) + 0.5*dH_N2(i)) ; enddo
      ! if ((.not.N2_region_set(i)) .and. (dH_N2(i) > 0.5*dH_subML)) then
      !    ! Use whatever stratification we can, measured over whatever distance is available?
      !    T_deeper(i) = tv%T(i,j,nz) ; S_deeper(i) = tv%S(i,j,nz)
      !    N2_region_set(i) = .true.
      ! endif
      call calculate_density(T_subML, S_subML, pRef_N2, rho_subML, tv%eqn_of_state, EOSdom)
      call calculate_density(T_deeper, S_deeper, pRef_N2, rho_deeper, tv%eqn_of_state, EOSdom)
      do i=is,ie ; if ((G%mask2dT(i,j)>0.5) .and. N2_region_set(i)) then
        subMLN2(i,j) =  gE_rho0 * (rho_deeper(i) - rho_subML(i)) / (GV%H_to_z * dH_N2(i))
      endif ; enddo
    endif
  enddo ! j-loop

  if (id_MLD > 0) call post_data(id_MLD, MLD, diagPtr)
  if (id_N2 > 0)  call post_data(id_N2, subMLN2, diagPtr)
  if (id_SQ > 0)  call post_data(id_SQ, MLD2, diagPtr)

end subroutine diagnoseMLDbyDensityDifference

!> Diagnose a mixed layer depth (MLD) determined by the depth a given energy value would mix.
!> This routine is appropriate in MOM_diabatic_driver due to its position within the time stepping.
subroutine diagnoseMLDbyEnergy(id_MLD, h, tv, G, GV, US, Mixing_Energy, diagPtr)
  ! Author: Brandon Reichl
  ! Date: October 2, 2020
  ! //
  ! *Note that gravity is assumed constant everywhere and divided out of all calculations.
  !
  ! This code has been written to step through the columns layer by layer, summing the PE
  ! change inferred by mixing the layer with all layers above.  When the change exceeds a
  ! threshold (determined by input array Mixing_Energy), the code needs to solve for how far
  ! into this layer the threshold PE change occurs (assuming constant density layers).
  ! This is expressed here via solving the function F(X) = 0 where:
  ! F(X) = 0.5 * ( Ca*X^3/(D1+X) + Cb*X^2/(D1+X) + Cc*X/(D1+X) + Dc/(D1+X)
  !                + Ca2*X^2 + Cb2*X + Cc2)
  ! where all coefficients are determined by the previous mixed layer depth, the
  ! density of the previous mixed layer, the present layer thickness, and the present
  ! layer density.  This equation is worked out by computing the total PE assuming constant
  ! density in the mixed layer as well as in the remaining part of the present layer that is
  ! not mixed.
  ! To solve for X in this equation a Newton's method iteration is employed, which
  ! converges extremely quickly (usually 1 guess) since this equation turns out being rather
  ! lienar for PE change with increasing X.
  ! Input parameters:
  integer, dimension(3),   intent(in) :: id_MLD      !< Energy output diag IDs
  type(ocean_grid_type),   intent(in) :: G           !< Grid type
  type(verticalGrid_type), intent(in) :: GV          !< ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US          !< A dimensional unit scaling type
  real, dimension(3),      intent(in) :: Mixing_Energy !< Energy values for up to 3 MLDs
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h           !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in) :: tv          !< Structure containing pointers to any
                                                     !! available thermodynamic fields.
  type(diag_ctrl),         pointer    :: diagPtr     !< Diagnostics structure

  ! Local variables
  real, dimension(SZI_(G), SZJ_(G),3) :: MLD     ! Diagnosed mixed layer depth [Z ~> m].
  real, dimension(SZK_(G)) :: Z_L, Z_U, dZ, Rho_c, pRef_MLD
  real, dimension(3) :: PE_threshold

  real :: ig, E_g
  real :: PE_Threshold_fraction, PE, PE_Mixed, PE_Mixed_TST
  real :: RhoDZ_ML, H_ML, RhoDZ_ML_TST, H_ML_TST
  real :: Rho_ML

  real :: R1, D1, R2, D2
  real :: Ca, Cb,D ,Cc, Cd, Ca2, Cb2, C, Cc2
  real :: Gx, Gpx, Hx, iHx, Hpx, Ix, Ipx, Fgx, Fpx, X, X2

  integer :: IT, iM
  integer :: i, j, is, ie, js, je, k, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  pRef_MLD(:) = 0.0
  mld(:,:,:) = 0.0
  PE_Threshold_fraction = 1.e-4 !Fixed threshold of 0.01%, could be runtime.

  do iM=1,3
    PE_threshold(iM) = Mixing_Energy(iM)/GV%g_earth
  enddo

  do j=js,je; do i=is,ie
    if (G%mask2dT(i,j) > 0.0) then

      call calculate_density(tv%T(i,j,:), tv%S(i,j,:), pRef_MLD, rho_c, 1, nz, &
                             tv%eqn_of_state, scale=US%kg_m3_to_R)

      do k=1,nz
        DZ(k) = h(i,j,k) * GV%H_to_Z
      enddo
      Z_U(1) = 0.0
      Z_L(1) = -DZ(1)
      do k=2,nz
        Z_U(k) = Z_L(k-1)
        Z_L(k) = Z_L(k-1)-DZ(k)
      enddo

      do iM=1,3

        ! Initialize these for each column-wise calculation
        PE = 0.0
        RhoDZ_ML = 0.0
        H_ML = 0.0
        RhoDZ_ML_TST = 0.0
        H_ML_TST = 0.0
        PE_Mixed = 0.0

        do k=1,nz

          ! This is the unmixed PE cumulative sum from top down
          PE = PE + 0.5 * rho_c(k) * (Z_U(k)**2 - Z_L(k)**2)

          ! This is the depth and integral of density
          H_ML_TST = H_ML + DZ(k)
          RhoDZ_ML_TST = RhoDZ_ML + rho_c(k) * DZ(k)

          ! The average density assuming all layers including this were mixed
          Rho_ML = RhoDZ_ML_TST/H_ML_TST

          ! The PE assuming all layers including this were mixed
          ! Note that 0. could be replaced with "Surface", which doesn't have to be 0
          ! but 0 is a good reference value.
          PE_Mixed_TST = 0.5 * Rho_ML * (0.**2 - (0. - H_ML_TST)**2)

          ! Check if we supplied enough energy to mix to this layer
          if (PE_Mixed_TST - PE <= PE_threshold(iM)) then
            H_ML = H_ML_TST
            RhoDZ_ML = RhoDZ_ML_TST

          else ! If not, we need to solve where the energy ran out
            ! This will be done with a Newton's method iteration:

            R1 = RhoDZ_ML / H_ML ! The density of the mixed layer (not including this layer)
            D1 = H_ML ! The thickness of the mixed layer (not including this layer)
            R2 = rho_c(k) ! The density of this layer
            D2 = DZ(k) ! The thickness of this layer

            ! This block could be used to calculate the function coefficients if
            ! we don't reference all values to a surface designated as z=0
            ! S = Surface
            ! Ca  = -(R2)
            ! Cb  = -( (R1*D1) + R2*(2.*D1-2.*S) )
            ! D   = D1**2. - 2.*D1*S
            ! Cc  = -( R1*D1*(2.*D1-2.*S) + R2*D )
            ! Cd  = -(R1*D1*D)
            ! Ca2 = R2
            ! Cb2 = R2*(2*D1-2*S)
            ! C   = S**2 + D2**2 + D1**2 - 2*D1*S - 2.*D2*S +2.*D1*D2
            ! Cc2 = R2*(D+S**2-C)
            !
            ! If the surface is S = 0, it simplifies to:
            Ca  = -R2
            Cb  = -(R1 * D1 + R2 * (2. * D1))
            D   = D1**2
            Cc  = -(R1 * D1 * (2. * D1) + (R2 * D))
            Cd  = -R1 * (D1 * D)
            Ca2 = R2
            Cb2 = R2 * (2. * D1)
            C   = D2**2 + D1**2 + 2. * (D1 * D2)
            Cc2 = R2 * (D - C)

            ! First guess for an iteration using Newton's method
            X = DZ(k) * 0.5

            IT=0
            do while(IT<10)!We can iterate up to 10 times
              ! We are trying to solve the function:
              ! F(x) = G(x)/H(x)+I(x)
              ! for where F(x) = PE+PE_threshold, or equivalently for where
              ! F(x) = G(x)/H(x)+I(x) - (PE+PE_threshold) = 0
              ! We also need the derivative of this function for the Newton's method iteration
              ! F'(x) = (G'(x)H(x)-G(x)H'(x))/H(x)^2 + I'(x)
              ! G and its derivative
              Gx = 0.5 * (Ca * (X*X*X) + Cb * X**2 + Cc * X + Cd)
              Gpx = 0.5 * (3. * (Ca * X**2) + 2. * (Cb * X) + Cc)
              ! H, its inverse, and its derivative
              Hx = D1 + X
              iHx = 1. / Hx
              Hpx = 1.
              ! I and its derivative
              Ix = 0.5 * (Ca2 * X**2 + Cb2 * X + Cc2)
              Ipx = 0.5 * (2. * Ca2 * X + Cb2)

              ! The Function and its derivative:
              PE_Mixed = Gx * iHx + Ix
              Fgx = PE_Mixed - (PE + PE_threshold(iM))
              Fpx = (Gpx * Hx - Hpx * Gx) * iHx**2 + Ipx

              ! Check if our solution is within the threshold bounds, if not update
              ! using Newton's method.  This appears to converge almost always in
              ! one step because the function is very close to linear in most applications.
              if (abs(Fgx) > PE_Threshold(iM) * PE_Threshold_fraction) then
                X2 = X - Fgx / Fpx
                IT = IT + 1
                if (X2 < 0. .or. X2 > DZ(k)) then
                  ! The iteration seems to be robust, but we need to do something *if*
                  ! things go wrong... How should we treat failed iteration?
                  ! Present solution: Stop trying to compute and just say we can't mix this layer.
                  X=0
                  exit
                else
                  X = X2
                endif
              else
                exit! Quit the iteration
              endif
            enddo
            H_ML = H_ML + X
            exit! Quit looping through the column
          endif
        enddo
        MLD(i,j,iM) = H_ML
      enddo
    else
      MLD(i,j,:) = 0.0
    endif
  enddo ; enddo

  if (id_MLD(1) > 0) call post_data(id_MLD(1), MLD(:,:,1), diagPtr)
  if (id_MLD(2) > 0) call post_data(id_MLD(2), MLD(:,:,2), diagPtr)
  if (id_MLD(3) > 0) call post_data(id_MLD(3), MLD(:,:,3), diagPtr)

end subroutine diagnoseMLDbyEnergy

!> Update the thickness, temperature, and salinity due to thermodynamic
!! boundary forcing (contained in fluxes type) applied to h, tv%T and tv%S,
!! and calculate the TKE implications of this heating.
subroutine applyBoundaryFluxesInOut(CS, G, GV, US, dt, fluxes, optics, nsw, h, tv, &
                                    aggregate_FW_forcing, evap_CFL_limit, &
                                    minimum_forcing_depth, cTKE, dSV_dT, dSV_dS, &
                                    SkinBuoyFlux )
  type(diabatic_aux_CS),   pointer       :: CS !< Control structure for diabatic_aux
  type(ocean_grid_type),   intent(in)    :: G  !< Grid structure
  type(verticalGrid_type), intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real,                    intent(in)    :: dt !< Time-step over which forcing is applied [T ~> s]
  type(forcing),           intent(inout) :: fluxes !< Surface fluxes container
  type(optics_type),       pointer       :: optics !< Optical properties container
  integer,                 intent(in)    :: nsw !< The number of frequency bands of penetrating
                                                !! shortwave radiation
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv !< Structure containing pointers to any
                                               !! available thermodynamic fields.
  logical,                 intent(in)    :: aggregate_FW_forcing !< If False, treat in/out fluxes separately.
  real,                    intent(in)    :: evap_CFL_limit !< The largest fraction of a layer that
                                               !! can be evaporated in one time-step [nondim].
  real,                    intent(in)    :: minimum_forcing_depth !< The smallest depth over which
                                               !! heat and freshwater fluxes is applied [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                 optional, intent(out)   :: cTKE !< Turbulent kinetic energy requirement to mix
                                               !! forcing through each layer [R Z3 T-2 ~> J m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                 optional, intent(out)   :: dSV_dT !< Partial derivative of specific volume with
                                               !! potential temperature [R-1 degC-1 ~> m3 kg-1 degC-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                 optional, intent(out)   :: dSV_dS !< Partial derivative of specific volume with
                                               !! salinity [R-1 ppt-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZI_(G),SZJ_(G)), &
                   optional, intent(out) :: SkinBuoyFlux !< Buoyancy flux at surface [Z2 T-3 ~> m2 s-3].

  ! Local variables
  integer, parameter :: maxGroundings = 5
  integer :: numberOfGroundings, iGround(maxGroundings), jGround(maxGroundings)
  real :: H_limit_fluxes
  real :: IforcingDepthScale
  real :: Idt        ! The inverse of the timestep [T-1 ~> s-1]
  real :: dThickness, dTemp, dSalt
  real :: fractionOfForcing, hOld, Ithickness
  real :: RivermixConst  ! A constant used in implementing river mixing [R Z2 T-1 ~> Pa s].

  real, dimension(SZI_(G)) :: &
    d_pres,       &  ! pressure change across a layer [R L2 T-2 ~> Pa]
    p_lay,        &  ! average pressure in a layer [R L2 T-2 ~> Pa]
    pres,         &  ! pressure at an interface [R L2 T-2 ~> Pa]
    netMassInOut, &  ! surface water fluxes [H ~> m or kg m-2] over time step
    netMassIn,    &  ! mass entering ocean surface [H ~> m or kg m-2] over a time step
    netMassOut,   &  ! mass leaving ocean surface [H ~> m or kg m-2] over a time step
    netHeat,      &  ! heat via surface fluxes excluding Pen_SW_bnd and netMassOut
                     ! [degC H ~> degC m or degC kg m-2]
    netSalt,      &  ! surface salt flux ( g(salt)/m2 for non-Bouss and ppt*H for Bouss )
                     ! [ppt H ~> ppt m or ppt kg m-2]
    nonpenSW,     &  ! non-downwelling SW, which is absorbed at ocean surface
                     ! [degC H ~> degC m or degC kg m-2]
    SurfPressure, &  ! Surface pressure (approximated as 0.0) [R L2 T-2 ~> Pa]
    dRhodT,       &  ! change in density per change in temperature [R degC-1 ~> kg m-3 degC-1]
    dRhodS,       &  ! change in density per change in salinity [R ppt-1 ~> kg m-3 ppt-1]
    netheat_rate, &  ! netheat but for dt=1 [degC H T-1 ~> degC m s-1 or degC kg m-2 s-1]
    netsalt_rate, &  ! netsalt but for dt=1 (e.g. returns a rate)
                     ! [ppt H T-1 ~> ppt m s-1 or ppt kg m-2 s-1]
    netMassInOut_rate! netmassinout but for dt=1 [H T-1 ~> m s-1 or kg m-2 s-1]
  real, dimension(SZI_(G), SZK_(G)) :: &
    h2d, &           ! A 2-d copy of the thicknesses [H ~> m or kg m-2]
    T2d, &           ! A 2-d copy of the layer temperatures [degC]
    pen_TKE_2d, &    ! The TKE required to homogenize the heating by shortwave radiation within
                     ! a layer [R Z3 T-2 ~> J m-2]
    dSV_dT_2d        ! The partial derivative of specific volume with temperature [R-1 degC-1 ~> m3 kg-1 degC-1]
  real, dimension(SZI_(G)) :: &
    netPen_rate      ! The surface penetrative shortwave heating rate summed over all bands
                     ! [degC H T-1 ~> degC m s-1 or degC kg m-2 s-1]
  real, dimension(max(nsw,1),SZI_(G)) :: &
    Pen_SW_bnd, &    ! The penetrative shortwave heating integrated over a timestep by band
                     ! [degC H ~> degC m or degC kg m-2]
    Pen_SW_bnd_rate  ! The penetrative shortwave heating rate by band
                     ! [degC H T-1 ~> degC m s-1 or degC kg m-2 s-1]
  real, dimension(max(nsw,1),SZI_(G),SZK_(G)) :: &
    opacityBand      ! The opacity (inverse of the exponential absorption length) of each frequency
                     ! band of shortwave radation in each layer [H-1 ~> m-1 or m2 kg-1]
  real, dimension(maxGroundings) :: hGrounding ! Thickness added by each grounding event [H ~> m or kg m-2]
  real    :: Temp_in, Salin_in
  real    :: g_Hconv2 ! A conversion factor for use in the TKE calculation
                      ! in units of [Z3 R2 T-2 H-2 ~> kg2 m-5 s-2 or m s-2].
  real    :: GoRho    ! g_Earth times a unit conversion factor divided by density
                      ! [Z T-2 R-1 ~> m4 s-2 kg-1]
  logical :: calculate_energetics
  logical :: calculate_buoyancy
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, is, ie, js, je, k, nz, n, nb
  character(len=45) :: mesg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  Idt = 1.0 / dt

  calculate_energetics = (present(cTKE) .and. present(dSV_dT) .and. present(dSV_dS))
  calculate_buoyancy = present(SkinBuoyFlux)
  if (calculate_buoyancy) SkinBuoyFlux(:,:) = 0.0
  if (present(cTKE)) cTKE(:,:,:) = 0.0
  g_Hconv2 = (US%L_to_Z**2*GV%g_Earth * GV%H_to_RZ) * GV%H_to_RZ
  EOSdom(:) = EOS_domain(G%HI)

  ! Only apply forcing if fluxes%sw is associated.
  if (.not.associated(fluxes%sw) .and. .not.calculate_energetics) return

  if (calculate_buoyancy) then
    SurfPressure(:) = 0.0
    GoRho = US%L_to_Z**2*GV%g_Earth / GV%Rho0
  endif

  ! H_limit_fluxes is used by extractFluxes1d to scale down fluxes if the total
  ! depth of the ocean is vanishing. It does not (yet) handle a value of zero.
  ! To accommodate vanishing upper layers, we need to allow for an instantaneous
  ! distribution of forcing over some finite vertical extent. The bulk mixed layer
  ! code handles this issue properly.
  H_limit_fluxes = max(GV%Angstrom_H, 1.E-30*GV%m_to_H)

  ! diagnostic to see if need to create mass to avoid grounding
  if (CS%id_createdH>0) CS%createdH(:,:) = 0.
  numberOfGroundings = 0

  !$OMP parallel do default(none) shared(is,ie,js,je,nz,h,tv,nsw,G,GV,US,optics,fluxes,    &
  !$OMP                                  H_limit_fluxes,numberOfGroundings,iGround,jGround,&
  !$OMP                                  nonPenSW,hGrounding,CS,Idt,aggregate_FW_forcing,  &
  !$OMP                                  minimum_forcing_depth,evap_CFL_limit,dt,EOSdom,   &
  !$OMP                                  calculate_buoyancy,netPen_rate,SkinBuoyFlux,GoRho, &
  !$OMP                                  calculate_energetics,dSV_dT,dSV_dS,cTKE,g_Hconv2) &
  !$OMP                          private(opacityBand,h2d,T2d,netMassInOut,netMassOut,      &
  !$OMP                                  netHeat,netSalt,Pen_SW_bnd,fractionOfForcing,     &
  !$OMP                                  IforcingDepthScale,                               &
  !$OMP                                  dThickness,dTemp,dSalt,hOld,Ithickness,           &
  !$OMP                                  netMassIn,pres,d_pres,p_lay,dSV_dT_2d,            &
  !$OMP                                  netmassinout_rate,netheat_rate,netsalt_rate,      &
  !$OMP                                  drhodt,drhods,pen_sw_bnd_rate,                    &
  !$OMP                                  pen_TKE_2d,Temp_in,Salin_in,RivermixConst)        &
  !$OMP                     firstprivate(SurfPressure)
  do j=js,je
  ! Work in vertical slices for efficiency

    ! Copy state into 2D-slice arrays
    do k=1,nz ; do i=is,ie
      h2d(i,k) = h(i,j,k)
      T2d(i,k) = tv%T(i,j,k)
    enddo ; enddo

    if (calculate_energetics) then
      ! The partial derivatives of specific volume with temperature and
      ! salinity need to be precalculated to avoid having heating of
      ! tiny layers give nonsensical values.
      if (associated(tv%p_surf)) then
        do i=is,ie ; pres(i) = tv%p_surf(i,j) ; enddo
      else
        do i=is,ie ; pres(i) = 0.0 ; enddo
      endif
      do k=1,nz
        do i=is,ie
          d_pres(i) = (GV%g_Earth * GV%H_to_RZ) * h2d(i,k)
          p_lay(i) = pres(i) + 0.5*d_pres(i)
          pres(i) = pres(i) + d_pres(i)
        enddo
        call calculate_specific_vol_derivs(T2d(:,k), tv%S(:,j,k), p_lay(:), &
                 dSV_dT(:,j,k), dSV_dS(:,j,k), tv%eqn_of_state, EOSdom)
        do i=is,ie ; dSV_dT_2d(i,k) = dSV_dT(i,j,k) ; enddo
      enddo
      pen_TKE_2d(:,:) = 0.0
    endif

    ! Nothing more is done on this j-slice if there is no buoyancy forcing.
    if (.not.associated(fluxes%sw)) cycle

    if (nsw>0) call extract_optics_slice(optics, j, G, GV, opacity=opacityBand, opacity_scale=(1.0/GV%m_to_H))

    ! The surface forcing is contained in the fluxes type.
    ! We aggregate the thermodynamic forcing for a time step into the following:
    ! netMassInOut = surface water fluxes [H ~> m or kg m-2] over time step
    !              = lprec + fprec + vprec + evap + lrunoff + frunoff
    !                note that lprec generally has sea ice melt/form included.
    ! netMassOut   = net mass leaving ocean surface [H ~> m or kg m-2] over a time step.
    !                netMassOut < 0 means mass leaves ocean.
    ! netHeat      = heat via surface fluxes [degC H ~> degC m or degC kg m-2], excluding the part
    !                contained in Pen_SW_bnd; and excluding heat_content of netMassOut < 0.
    ! netSalt      = surface salt fluxes [ppt H ~> dppt m or gSalt m-2]
    ! Pen_SW_bnd   = components to penetrative shortwave radiation split according to bands.
    !                This field provides that portion of SW from atmosphere that in fact
    !                enters to the ocean and participates in pentrative SW heating.
    ! nonpenSW     = non-downwelling SW flux, which is absorbed in ocean surface
    !                (in tandem w/ LW,SENS,LAT); saved only for diagnostic purposes.

    !----------------------------------------------------------------------------------------
    !BGR-June 26, 2017{
    !Temporary action to preserve answers while fixing a bug.
    ! To fix a bug in a diagnostic calculation, applyboundaryfluxesinout now returns
    !  the surface buoyancy flux. Previously, extractbuoyancyflux2d was called, meaning
    !  a second call to extractfluxes1d (causing the diagnostic net_heat to be incorrect).
    !  Note that this call to extract buoyancyflux2d was AFTER applyboundaryfluxesinout,
    !  which means it used the T/S fields after this routine.  Therefore, the surface
    !  buoyancy flux is computed here at the very end of this routine for legacy reasons.
    !  A few specific notes follow:
    !     1) The old method did not included river/calving contributions to heat flux.  This
    !        is kept consistent here via commenting code in the present extractFluxes1d <_rate>
    !        outputs, but we may reconsider this approach.
    !     2) The old method computed the buoyancy flux rate directly (by setting dt=1), instead
    !        of computing the integrated value (and dividing by dt). Hence the required
    !        additional outputs from extractFluxes1d.
    !          *** This is because: A*dt/dt =/=  A due to round off.
    !     3) The old method computed buoyancy flux after this routine, meaning the returned
    !        surface fluxes (from extractfluxes1d) must be recorded for use later in the code.
    !        We could (and maybe should) move that loop up to before the surface fluxes are
    !        applied, but this will change answers.
    !     For all these reasons we compute additional values of <_rate> which are preserved
    !     for the buoyancy flux calculation and reproduce the old answers.
    !   In the future this needs more detailed investigation to make sure everything is
    !   consistent and correct. These details shouldnt significantly effect climate,
    !   but do change answers.
    !-----------------------------------------------------------------------------------------
    if (calculate_buoyancy) then
      call extractFluxes1d(G, GV, US, fluxes, optics, nsw, j, dt,          &
                  H_limit_fluxes, CS%use_river_heat_content, CS%use_calving_heat_content, &
                  h2d, T2d, netMassInOut, netMassOut, netHeat, netSalt,                   &
                  Pen_SW_bnd, tv, aggregate_FW_forcing, nonpenSW=nonpenSW,                &
                  net_Heat_rate=netheat_rate, net_salt_rate=netsalt_rate,                 &
                  netmassinout_rate=netmassinout_rate, pen_sw_bnd_rate=pen_sw_bnd_rate)
    else
      call extractFluxes1d(G, GV, US, fluxes, optics, nsw, j, dt,          &
                  H_limit_fluxes, CS%use_river_heat_content, CS%use_calving_heat_content, &
                  h2d, T2d, netMassInOut, netMassOut, netHeat, netSalt,                   &
                  Pen_SW_bnd, tv, aggregate_FW_forcing, nonpenSW=nonpenSW)
   endif
    ! ea is for passive tracers
    do i=is,ie
    !  ea(i,j,1) = netMassInOut(i)
      if (aggregate_FW_forcing) then
        netMassOut(i) = netMassInOut(i)
        netMassIn(i) = 0.
      else
        netMassIn(i) = netMassInOut(i) - netMassOut(i)
      endif
      if (G%mask2dT(i,j)>0.0) then
        fluxes%netMassOut(i,j) = netMassOut(i)
        fluxes%netMassIn(i,j) = netMassIn(i)
      else
        fluxes%netMassOut(i,j) = 0.0
        fluxes%netMassIn(i,j) = 0.0
      endif
    enddo

    ! Apply the surface boundary fluxes in three steps:
    ! A/ update mass, temp, and salinity due to all terms except mass leaving
    !    ocean (and corresponding outward heat content), and ignoring penetrative SW.
    ! B/ update mass, salt, temp from mass leaving ocean.
    ! C/ update temp due to penetrative SW
    do i=is,ie
      if (G%mask2dT(i,j)>0.) then

        ! A/ Update mass, temp, and salinity due to incoming mass flux.
        do k=1,1

          ! Change in state due to forcing
          dThickness = netMassIn(i) ! Since we are adding mass, we can use all of it
          dTemp = 0.
          dSalt = 0.

          ! Update the forcing by the part to be consumed within the present k-layer.
          ! If fractionOfForcing = 1, then updated netMassIn, netHeat, and netSalt vanish.
          netMassIn(i) = netMassIn(i) - dThickness
          ! This line accounts for the temperature of the mass exchange
          Temp_in = T2d(i,k)
          Salin_in = 0.0
          dTemp = dTemp + dThickness*Temp_in

          ! Diagnostics of heat content associated with mass fluxes
          if (associated(fluxes%heat_content_massin))                             &
            fluxes%heat_content_massin(i,j) = fluxes%heat_content_massin(i,j) +   &
                         T2d(i,k) * max(0.,dThickness) * GV%H_to_RZ * fluxes%C_p * Idt
          if (associated(fluxes%heat_content_massout))                            &
            fluxes%heat_content_massout(i,j) = fluxes%heat_content_massout(i,j) + &
                         T2d(i,k) * min(0.,dThickness) * GV%H_to_RZ * fluxes%C_p * Idt
          if (associated(tv%TempxPmE)) tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + &
                         T2d(i,k) * dThickness * GV%H_to_RZ

          ! Determine the energetics of river mixing before updating the state.
          if (calculate_energetics .and. associated(fluxes%lrunoff) .and. CS%do_rivermix) then
            ! Here we add an additional source of TKE to the mixed layer where river
            ! is present to simulate unresolved estuaries. The TKE input, TKE_river in
            ! [Z3 T-3 ~> m3 s-3], is diagnosed as follows:
            !   TKE_river = 0.5*rivermix_depth*g*(1/rho)*drho_ds*
            !               River*(Samb - Sriver) = CS%mstar*U_star^3
            ! where River is in units of [Z T-1 ~> m s-1].
            ! Samb = Ambient salinity at the mouth of the estuary
            ! rivermix_depth =  The prescribed depth over which to mix river inflow
            ! drho_ds = The gradient of density wrt salt at the ambient surface salinity.
            ! Sriver = 0 (i.e. rivers are assumed to be pure freshwater)
            if (GV%Boussinesq) then
              RivermixConst = -0.5*(CS%rivermix_depth*dt) * ( US%L_to_Z**2*GV%g_Earth ) * GV%Rho0
            else
              RivermixConst = -0.5*(CS%rivermix_depth*dt) * GV%Rho0 * ( US%L_to_Z**2*GV%g_Earth )
            endif
            cTKE(i,j,k) = cTKE(i,j,k) + max(0.0, RivermixConst*dSV_dS(i,j,1) * &
                            (fluxes%lrunoff(i,j) + fluxes%frunoff(i,j)) * tv%S(i,j,1))
          endif

          ! Update state
          hOld     = h2d(i,k)               ! Keep original thickness in hand
          h2d(i,k) = h2d(i,k) + dThickness  ! New thickness
          if (h2d(i,k) > 0.0) then
            if (calculate_energetics .and. (dThickness > 0.)) then
              ! Calculate the energy required to mix the newly added water over
              ! the topmost grid cell.
              cTKE(i,j,k) = cTKE(i,j,k) + 0.5*g_Hconv2*(hOld*dThickness) * &
                 ((T2d(i,k) - Temp_in) * dSV_dT(i,j,k) + (tv%S(i,j,k) - Salin_in) * dSV_dS(i,j,k))
            endif
            Ithickness  = 1.0/h2d(i,k)      ! Inverse new thickness
            ! The "if"s below avoid changing T/S by roundoff unnecessarily
            if (dThickness /= 0. .or. dTemp /= 0.) T2d(i,k)    = (hOld*T2d(i,k)    + dTemp)*Ithickness
            if (dThickness /= 0. .or. dSalt /= 0.) tv%S(i,j,k) = (hOld*tv%S(i,j,k) + dSalt)*Ithickness

          endif

        enddo ! k=1,1

        ! B/ Update mass, salt, temp from mass leaving ocean and other fluxes of heat and salt.
        do k=1,nz

          ! Place forcing into this layer if this layer has nontrivial thickness.
          ! For layers thin relative to 1/IforcingDepthScale, then distribute
          ! forcing into deeper layers.
          IforcingDepthScale = 1. / max(GV%H_subroundoff, minimum_forcing_depth - netMassOut(i) )
          ! fractionOfForcing = 1.0, unless h2d is less than IforcingDepthScale.
          fractionOfForcing = min(1.0, h2d(i,k)*IforcingDepthScale)

          ! In the case with (-1)*netMassOut*fractionOfForcing greater than cfl*h, we
          ! limit the forcing applied to this cell, leaving the remaining forcing to
          ! be distributed downwards.
          if (-fractionOfForcing*netMassOut(i) > evap_CFL_limit*h2d(i,k)) then
            fractionOfForcing = -evap_CFL_limit*h2d(i,k)/netMassOut(i)
          endif

          ! Change in state due to forcing

          dThickness = max( fractionOfForcing*netMassOut(i), -h2d(i,k) )
          dTemp      = fractionOfForcing*netHeat(i)
          !   ### The 0.9999 here should become a run-time parameter?
          dSalt = max( fractionOfForcing*netSalt(i), -0.9999*h2d(i,k)*tv%S(i,j,k))

          ! Update the forcing by the part to be consumed within the present k-layer.
          ! If fractionOfForcing = 1, then new netMassOut vanishes.
          netMassOut(i) = netMassOut(i) - dThickness
          netHeat(i) = netHeat(i) - dTemp
          netSalt(i) = netSalt(i) - dSalt

          ! This line accounts for the temperature of the mass exchange
          dTemp = dTemp + dThickness*T2d(i,k)

          ! Diagnostics of heat content associated with mass fluxes
          if (associated(fluxes%heat_content_massin)) &
            fluxes%heat_content_massin(i,j) = fluxes%heat_content_massin(i,j) + &
                         T2d(i,k) * max(0.,dThickness) * GV%H_to_RZ * fluxes%C_p * Idt
          if (associated(fluxes%heat_content_massout)) &
            fluxes%heat_content_massout(i,j) = fluxes%heat_content_massout(i,j) + &
                         T2d(i,k) * min(0.,dThickness) * GV%H_to_RZ * fluxes%C_p * Idt
          if (associated(tv%TempxPmE)) tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + &
                         T2d(i,k) * dThickness * GV%H_to_RZ

          ! Update state by the appropriate increment.
          hOld     = h2d(i,k)               ! Keep original thickness in hand
          h2d(i,k) = h2d(i,k) + dThickness  ! New thickness

          if (h2d(i,k) > 0.) then
            if (calculate_energetics) then
              ! Calculate the energy required to mix the newly added water over the topmost grid
              ! cell, assuming that the fluxes of heat and salt and rejected brine are initially
              ! applied in vanishingly thin layers at the top of the layer before being mixed
              ! throughout the layer.  Note that dThickness is always <= 0 here, and that
              ! negative cTKE is a deficit that will need to be filled later.
              cTKE(i,j,k) = cTKE(i,j,k) - (0.5*h2d(i,k)*g_Hconv2) * &
                            ((dTemp - dthickness*T2d(i,k)) * dSV_dT(i,j,k) + &
                             (dSalt - dthickness*tv%S(i,j,k)) * dSV_dS(i,j,k))
            endif
            Ithickness  = 1.0/h2d(i,k) ! Inverse of new thickness
            T2d(i,k)    = (hOld*T2d(i,k) + dTemp)*Ithickness
            tv%S(i,j,k) = (hOld*tv%S(i,j,k) + dSalt)*Ithickness
          elseif (h2d(i,k) < 0.0) then ! h2d==0 is a special limit that needs no extra handling
            call forcing_SinglePointPrint(fluxes,G,i,j,'applyBoundaryFluxesInOut (h<0)')
            write(0,*) 'applyBoundaryFluxesInOut(): lon,lat=',G%geoLonT(i,j),G%geoLatT(i,j)
            write(0,*) 'applyBoundaryFluxesInOut(): netT,netS,netH=',netHeat(i),netSalt(i),netMassInOut(i)
            write(0,*) 'applyBoundaryFluxesInOut(): dT,dS,dH=',dTemp,dSalt,dThickness
            write(0,*) 'applyBoundaryFluxesInOut(): h(n),h(n+1),k=',hOld,h2d(i,k),k
            call MOM_error(FATAL, "MOM_diabatic_driver.F90, applyBoundaryFluxesInOut(): "//&
                           "Complete mass loss in column!")
          endif

        enddo ! k

      ! Check if trying to apply fluxes over land points
      elseif ((abs(netHeat(i))+abs(netSalt(i))+abs(netMassIn(i))+abs(netMassOut(i)))>0.) then

        if (.not. CS%ignore_fluxes_over_land) then
           call forcing_SinglePointPrint(fluxes,G,i,j,'applyBoundaryFluxesInOut (land)')
           write(0,*) 'applyBoundaryFluxesInOut(): lon,lat=',G%geoLonT(i,j),G%geoLatT(i,j)
           write(0,*) 'applyBoundaryFluxesInOut(): netHeat,netSalt,netMassIn,netMassOut=',&
                   netHeat(i),netSalt(i),netMassIn(i),netMassOut(i)

           call MOM_error(FATAL, "MOM_diabatic_driver.F90, applyBoundaryFluxesInOut(): "//&
                                 "Mass loss over land?")
        endif

      endif

      ! If anything remains after the k-loop, then we have grounded out, which is a problem.
      if (netMassIn(i)+netMassOut(i) /= 0.0) then
!$OMP critical
        numberOfGroundings = numberOfGroundings +1
        if (numberOfGroundings<=maxGroundings) then
          iGround(numberOfGroundings) = i ! Record i,j location of event for
          jGround(numberOfGroundings) = j ! warning message
          hGrounding(numberOfGroundings) = netMassIn(i)+netMassOut(i)
        endif
!$OMP end critical
        if (CS%id_createdH>0) CS%createdH(i,j) = CS%createdH(i,j) - (netMassIn(i)+netMassOut(i))/dt
      endif

    enddo ! i

    ! Step C/ in the application of fluxes
    ! Heat by the convergence of penetrating SW.
    ! SW penetrative heating uses the updated thickness from above.

    ! Save temperature before increment with SW heating
    ! and initialize CS%penSWflux_diag to zero.
    if (CS%id_penSW_diag > 0 .or. CS%id_penSWflux_diag > 0) then
      do k=1,nz ; do i=is,ie
        CS%penSW_diag(i,j,k)     = T2d(i,k)
        CS%penSWflux_diag(i,j,k) = 0.0
      enddo ; enddo
      k=nz+1 ; do i=is,ie
        CS%penSWflux_diag(i,j,k) = 0.0
      enddo
    endif

    if (calculate_energetics) then
      call absorbRemainingSW(G, GV, US, h2d, opacityBand, nsw, optics, j, dt, H_limit_fluxes, &
                             .false., .true., T2d, Pen_SW_bnd, TKE=pen_TKE_2d, dSV_dT=dSV_dT_2d)
      k = 1 ! For setting break-points.
      do k=1,nz ; do i=is,ie
        cTKE(i,j,k) = cTKE(i,j,k) + pen_TKE_2d(i,k)
      enddo ; enddo
    else
      call absorbRemainingSW(G, GV, US, h2d, opacityBand, nsw, optics, j, dt, H_limit_fluxes, &
                             .false., .true., T2d, Pen_SW_bnd)
    endif


    ! Step D/ copy updated thickness and temperature
    ! 2d slice now back into model state.
    do k=1,nz ; do i=is,ie
      h(i,j,k)    = h2d(i,k)
      tv%T(i,j,k) = T2d(i,k)
    enddo ; enddo

    ! Diagnose heating [Q R Z T-1 ~> W m-2] applied to a grid cell from SW penetration
    ! Also diagnose the penetrative SW heat flux at base of layer.
    if (CS%id_penSW_diag > 0 .or. CS%id_penSWflux_diag > 0) then

      ! convergence of SW into a layer
      do k=1,nz ; do i=is,ie
        CS%penSW_diag(i,j,k) = (T2d(i,k)-CS%penSW_diag(i,j,k))*h(i,j,k) * Idt * tv%C_p * GV%H_to_RZ
      enddo ; enddo

      ! Perform a cumulative sum upwards from bottom to
      ! diagnose penetrative SW flux at base of tracer cell.
      ! CS%penSWflux_diag(i,j,k=1)    is penetrative shortwave at top of ocean.
      ! CS%penSWflux_diag(i,j,k=kbot+1) is zero, since assume no SW penetrates rock.
      ! CS%penSWflux_diag = rsdo  and CS%penSW_diag = rsdoabsorb
      ! rsdoabsorb(k) = rsdo(k) - rsdo(k+1), so that rsdo(k) = rsdo(k+1) + rsdoabsorb(k)
      if (CS%id_penSWflux_diag > 0) then
        do k=nz,1,-1 ; do i=is,ie
          CS%penSWflux_diag(i,j,k) = CS%penSW_diag(i,j,k) + CS%penSWflux_diag(i,j,k+1)
        enddo ; enddo
      endif

    endif

    ! Fill CS%nonpenSW_diag
    if (CS%id_nonpenSW_diag > 0) then
      do i=is,ie
        CS%nonpenSW_diag(i,j) = nonpenSW(i) * Idt * tv%C_p * GV%H_to_RZ
      enddo
    endif

    ! BGR: Get buoyancy flux to return for ePBL
    !  We want the rate, so we use the rate values returned from extractfluxes1d.
    !  Note that the *dt values could be divided by dt here, but
    !  1) Answers will change due to round-off
    !  2) Be sure to save their values BEFORE fluxes are used.
    if (Calculate_Buoyancy) then
      drhodt(:) = 0.0
      drhods(:) = 0.0
      netPen_rate(:) = 0.0
      ! Sum over bands and attenuate as a function of depth.
      ! netPen_rate is the netSW as a function of depth, but only the surface value is used here,
      ! in which case the values of dt, h, optics and H_limit_fluxes are irrelevant.  Consider
      ! writing a shorter and simpler variant to handle this very limited case.
      ! call sumSWoverBands(G, GV, US, h2d(:,:), optics_nbands(optics), optics, j, dt, &
      !                     H_limit_fluxes, .true., pen_SW_bnd_rate, netPen)
      do i=is,ie ; do nb=1,nsw ; netPen_rate(i) = netPen_rate(i) + pen_SW_bnd_rate(nb,i) ; enddo ; enddo

      ! Density derivatives
      if (associated(tv%p_surf)) then ; do i=is,ie ; SurfPressure(i) = tv%p_surf(i,j) ; enddo ; endif
      call calculate_density_derivs(T2d(:,1), tv%S(:,j,1), SurfPressure, dRhodT, dRhodS, &
                                    tv%eqn_of_state, EOSdom)
      ! 1. Adjust netSalt to reflect dilution effect of FW flux
      ! 2. Add in the SW heating for purposes of calculating the net
      ! surface buoyancy flux affecting the top layer.
      ! 3. Convert to a buoyancy flux, excluding penetrating SW heating
      !    BGR-Jul 5, 2017: The contribution of SW heating here needs investigated for ePBL.
      do i=is,ie
        SkinBuoyFlux(i,j) = - GoRho * GV%H_to_Z * &
            (dRhodS(i) * (netSalt_rate(i) - tv%S(i,j,1)*netMassInOut_rate(i)) + &
             dRhodT(i) * ( netHeat_rate(i) + netPen_rate(i)) ) ! [Z2 T-3 ~> m2 s-3]
      enddo
    endif

  enddo ! j-loop finish

  ! Post the diagnostics
  if (CS%id_createdH       > 0) call post_data(CS%id_createdH      , CS%createdH      , CS%diag)
  if (CS%id_penSW_diag     > 0) call post_data(CS%id_penSW_diag    , CS%penSW_diag    , CS%diag)
  if (CS%id_penSWflux_diag > 0) call post_data(CS%id_penSWflux_diag, CS%penSWflux_diag, CS%diag)
  if (CS%id_nonpenSW_diag  > 0) call post_data(CS%id_nonpenSW_diag , CS%nonpenSW_diag , CS%diag)

! The following check will be ignored if ignore_fluxes_over_land = true
  if (numberOfGroundings>0 .and. .not. CS%ignore_fluxes_over_land) then
    do i = 1, min(numberOfGroundings, maxGroundings)
      call forcing_SinglePointPrint(fluxes,G,iGround(i),jGround(i),'applyBoundaryFluxesInOut (grounding)')
      write(mesg(1:45),'(3es15.3)') G%geoLonT( iGround(i), jGround(i) ), &
                             G%geoLatT( iGround(i), jGround(i)), hGrounding(i)*GV%H_to_m
      call MOM_error(WARNING, "MOM_diabatic_driver.F90, applyBoundaryFluxesInOut(): "//&
                              "Mass created. x,y,dh= "//trim(mesg), all_print=.true.)
    enddo

    if (numberOfGroundings - maxGroundings > 0) then
      write(mesg, '(i4)') numberOfGroundings - maxGroundings
      call MOM_error(WARNING, "MOM_diabatic_driver:F90, applyBoundaryFluxesInOut(): "//&
                              trim(mesg) // " groundings remaining")
    endif
  endif

end subroutine applyBoundaryFluxesInOut

!> This subroutine initializes the parameters and control structure of the diabatic_aux module.
subroutine diabatic_aux_init(Time, G, GV, US, param_file, diag, CS, useALEalgorithm, use_ePBL)
  type(time_type), target, intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target, intent(inout) :: diag !< A structure used to regulate diagnostic output
  type(diabatic_aux_CS),   pointer       :: CS   !< A pointer to the control structure for the
                                                 !! diabatic_aux module, which is initialized here.
  logical,                 intent(in)    :: useALEalgorithm !< If true, use the ALE algorithm rather
                                                 !! than layered mode.
  logical,                 intent(in)    :: use_ePBL !< If true, use the implicit energetics planetary
                                                 !! boundary layer scheme to determine the diffusivity
                                                 !! in the surface boundary layer.

! This "include" declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl  = "MOM_diabatic_aux" ! This module's name.
  character(len=48)  :: thickness_units
  character(len=200) :: inputdir   ! The directory where NetCDF input files
  character(len=240) :: chl_filename ! A file from which chl_a concentrations are to be read.
  character(len=128) :: chl_file ! Data containing chl_a concentrations. Used
                                 ! when var_pen_sw is defined and reading from file.
  character(len=32)  :: chl_varname ! Name of chl_a variable in chl_file.
  logical :: use_temperature     ! True if thermodynamics are enabled.
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz, nbands
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "diabatic_aux_init called with an "// &
                            "associated control structure.")
    return
  else
    allocate(CS)
   endif

  CS%diag => diag
  CS%Time => Time

! Set default, read and log parameters
  call log_version(param_file, mdl, version, &
                   "The following parameters are used for auxiliary diabatic processes.")

  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", use_temperature, &
                 "If true, temperature and salinity are used as state "//&
                 "variables.", default=.true.)

  call get_param(param_file, mdl, "RECLAIM_FRAZIL", CS%reclaim_frazil, &
                 "If true, try to use any frazil heat deficit to cool any "//&
                 "overlying layers down to the freezing point, thereby "//&
                 "avoiding the creation of thin ice when the SST is above "//&
                 "the freezing point.", default=.true.)
  call get_param(param_file, mdl, "PRESSURE_DEPENDENT_FRAZIL", &
                                CS%pressure_dependent_frazil, &
                 "If true, use a pressure dependent freezing temperature "//&
                 "when making frazil. The default is false, which will be "//&
                 "faster but is inappropriate with ice-shelf cavities.", &
                 default=.false.)

  if (use_ePBL) then
    call get_param(param_file, mdl, "IGNORE_FLUXES_OVER_LAND", CS%ignore_fluxes_over_land,&
         "If true, the model does not check if fluxes are being applied "//&
         "over land points. This is needed when the ocean is coupled "//&
         "with ice shelves and sea ice, since the sea ice mask needs to "//&
         "be different than the ocean mask to avoid sea ice formation "//&
         "under ice shelves. This flag only works when use_ePBL = True.", default=.false.)
    call get_param(param_file, mdl, "DO_RIVERMIX", CS%do_rivermix, &
                 "If true, apply additional mixing wherever there is "//&
                 "runoff, so that it is mixed down to RIVERMIX_DEPTH "//&
                 "if the ocean is that deep.", default=.false.)
    if (CS%do_rivermix) &
      call get_param(param_file, mdl, "RIVERMIX_DEPTH", CS%rivermix_depth, &
                 "The depth to which rivers are mixed if DO_RIVERMIX is "//&
                 "defined.", units="m", default=0.0, scale=US%m_to_Z)
  else
    CS%do_rivermix = .false. ; CS%rivermix_depth = 0.0 ; CS%ignore_fluxes_over_land = .false.
  endif

  if (GV%nkml == 0) then
    call get_param(param_file, mdl, "USE_RIVER_HEAT_CONTENT", CS%use_river_heat_content, &
                   "If true, use the fluxes%runoff_Hflx field to set the "//&
                   "heat carried by runoff, instead of using SST*CP*liq_runoff.", &
                   default=.false.)
    call get_param(param_file, mdl, "USE_CALVING_HEAT_CONTENT", CS%use_calving_heat_content, &
                   "If true, use the fluxes%calving_Hflx field to set the "//&
                   "heat carried by runoff, instead of using SST*CP*froz_runoff.", &
                   default=.false.)
  else
    CS%use_river_heat_content = .false.
    CS%use_calving_heat_content = .false.
  endif

  if (useALEalgorithm) then
    CS%id_createdH = register_diag_field('ocean_model',"created_H",diag%axesT1, &
        Time, "The volume flux added to stop the ocean from drying out and becoming negative in depth", &
        "m s-1", conversion=GV%H_to_m*US%s_to_T)
    if (CS%id_createdH>0) allocate(CS%createdH(isd:ied,jsd:jed))

    ! diagnostic for heating of a grid cell from convergence of SW heat into the cell
    CS%id_penSW_diag = register_diag_field('ocean_model', 'rsdoabsorb',                     &
          diag%axesTL, Time, 'Convergence of Penetrative Shortwave Flux in Sea Water Layer',&
          'W m-2', conversion=US%QRZ_T_to_W_m2, &
          standard_name='net_rate_of_absorption_of_shortwave_energy_in_ocean_layer', v_extensive=.true.)

    ! diagnostic for penetrative SW heat flux at top interface of tracer cell (nz+1 interfaces)
    ! k=1 gives penetrative SW at surface; SW(k=nz+1)=0 (no penetration through rock).
    CS%id_penSWflux_diag = register_diag_field('ocean_model', 'rsdo',                               &
          diag%axesTi, Time, 'Downwelling Shortwave Flux in Sea Water at Grid Cell Upper Interface',&
          'W m-2', conversion=US%QRZ_T_to_W_m2, standard_name='downwelling_shortwave_flux_in_sea_water')

    ! need both arrays for the SW diagnostics (one for flux, one for convergence)
    if (CS%id_penSW_diag>0 .or. CS%id_penSWflux_diag>0) then
      allocate(CS%penSW_diag(isd:ied,jsd:jed,nz)) ; CS%penSW_diag(:,:,:) = 0.0
      allocate(CS%penSWflux_diag(isd:ied,jsd:jed,nz+1)) ; CS%penSWflux_diag(:,:,:) = 0.0
    endif

    ! diagnostic for non-downwelling SW radiation (i.e., SW absorbed at ocean surface)
    CS%id_nonpenSW_diag = register_diag_field('ocean_model', 'nonpenSW',                       &
          diag%axesT1, Time,                                                                   &
          'Non-downwelling SW radiation (i.e., SW absorbed in ocean surface with LW,SENS,LAT)',&
          'W m-2', conversion=US%QRZ_T_to_W_m2, &
          standard_name='nondownwelling_shortwave_flux_in_sea_water')
    if (CS%id_nonpenSW_diag > 0) then
      allocate(CS%nonpenSW_diag(isd:ied,jsd:jed)) ; CS%nonpenSW_diag(:,:) = 0.0
    endif
  endif

  if (use_temperature) then
    call get_param(param_file, mdl, "VAR_PEN_SW", CS%var_pen_sw, &
                   "If true, use one of the CHL_A schemes specified by "//&
                   "OPACITY_SCHEME to determine the e-folding depth of "//&
                   "incoming short wave radiation.", default=.false.)
    if (CS%var_pen_sw) then

      call get_param(param_file, mdl, "CHL_FROM_FILE", CS%chl_from_file, &
                   "If true, chl_a is read from a file.", default=.true.)
      if (CS%chl_from_file) then
        call time_interp_external_init()

        call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
        call get_param(param_file, mdl, "CHL_FILE", chl_file, &
                   "CHL_FILE is the file containing chl_a concentrations in "//&
                   "the variable CHL_A. It is used when VAR_PEN_SW and "//&
                   "CHL_FROM_FILE are true.", fail_if_missing=.true.)
        chl_filename = trim(slasher(inputdir))//trim(chl_file)
        call log_param(param_file, mdl, "INPUTDIR/CHL_FILE", chl_filename)
        call get_param(param_file, mdl, "CHL_VARNAME", chl_varname, &
                   "Name of CHL_A variable in CHL_FILE.", default='CHL_A')
        CS%sbc_chl = init_external_field(chl_filename, trim(chl_varname), domain=G%Domain%mpp_domain)
      endif

      CS%id_chl = register_diag_field('ocean_model', 'Chl_opac', diag%axesT1, Time, &
          'Surface chlorophyll A concentration used to find opacity', 'mg m-3')
    endif
  endif

  id_clock_uv_at_h = cpu_clock_id('(Ocean find_uv_at_h)', grain=CLOCK_ROUTINE)
  id_clock_frazil  = cpu_clock_id('(Ocean frazil)', grain=CLOCK_ROUTINE)

end subroutine diabatic_aux_init

!> This subroutine initializes the control structure and any related memory
!! for the diabatic_aux module.
subroutine diabatic_aux_end(CS)
  type(diabatic_aux_CS), pointer :: CS !< The control structure returned by a previous
                                       !! call to diabatic_aux_init; it is deallocated here.

  if (.not.associated(CS)) return

  if (CS%id_createdH       >0) deallocate(CS%createdH)
  if (CS%id_penSW_diag     >0) deallocate(CS%penSW_diag)
  if (CS%id_penSWflux_diag >0) deallocate(CS%penSWflux_diag)
  if (CS%id_nonpenSW_diag  >0) deallocate(CS%nonpenSW_diag)

  if (associated(CS)) deallocate(CS)

end subroutine diabatic_aux_end

!> \namespace mom_diabatic_aux
!!
!!    This module contains the subroutines that, along with the
!!  subroutines that it calls, implements diapycnal mass and momentum
!!  fluxes and a bulk mixed layer.  The diapycnal diffusion can be
!!  used without the bulk mixed layer.
!!
!!    diabatic first determines the (diffusive) diapycnal mass fluxes
!!  based on the convergence of the buoyancy fluxes within each layer.
!!  The dual-stream entrainment scheme of MacDougall and Dewar (JPO,
!!  1997) is used for combined diapycnal advection and diffusion,
!!  calculated implicitly and potentially with the Richardson number
!!  dependent mixing, as described by Hallberg (MWR, 2000). Diapycnal
!!  advection is fundamentally the residual of diapycnal diffusion,
!!  so the fully implicit upwind differencing scheme that is used is
!!  entirely appropriate.  The downward buoyancy flux in each layer
!!  is determined from an implicit calculation based on the previously
!!  calculated flux of the layer above and an estimated flux in the
!!  layer below.  This flux is subject to the following conditions:
!!  (1) the flux in the top and bottom layers are set by the boundary
!!  conditions, and (2) no layer may be driven below an Angstrom thick-
!!  ness.  If there is a bulk mixed layer, the buffer layer is treat-
!!  ed as a fixed density layer with vanishingly small diffusivity.
!!
!!    diabatic takes 5 arguments:  the two velocities (u and v), the
!!  thicknesses (h), a structure containing the forcing fields, and
!!  the length of time over which to act (dt).  The velocities and
!!  thickness are taken as inputs and modified within the subroutine.
!!  There is no limit on the time step.

end module MOM_diabatic_aux
