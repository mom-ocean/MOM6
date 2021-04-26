!> Calculates various values related to the bottom boundary layer, such as the viscosity and
!! thickness of the BBL (set_viscous_BBL).
module MOM_set_visc

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging, only : uvchksum, hchksum
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_domains, only : pass_var, CORNER
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, mech_forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_io, only : slasher, MOM_read_data
use MOM_kappa_shear, only : kappa_shear_is_used, kappa_shear_at_vertex
use MOM_cvmix_shear, only : cvmix_shear_is_used
use MOM_cvmix_conv,  only : cvmix_conv_is_used
use MOM_CVMix_ddiff, only : CVMix_ddiff_is_used
use MOM_restart, only : register_restart_field, query_initialized, MOM_restart_CS
use MOM_restart, only : register_restart_field_as_obsolete
use MOM_safe_alloc, only : safe_alloc_ptr, safe_alloc_alloc
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs, vertvisc_type
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, OBC_DIRECTION_E
use MOM_open_boundary, only : OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_open_boundary, only : OBC_segment_type
implicit none ; private

#include <MOM_memory.h>

public set_viscous_BBL, set_viscous_ML, set_visc_init, set_visc_end
public set_visc_register_restarts

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for MOM_set_visc
type, public :: set_visc_CS ; private
  real    :: Hbbl           !< The static bottom boundary layer thickness [H ~> m or kg m-2].
                            !! Runtime parameter `HBBL`.
  real    :: cdrag          !< The quadratic drag coefficient.
                            !! Runtime parameter `CDRAG`.
  real    :: c_Smag         !< The Laplacian Smagorinsky coefficient for
                            !! calculating the drag in channels.
  real    :: drag_bg_vel    !< An assumed unresolved background velocity for
                            !! calculating the bottom drag [L T-1 ~> m s-1].
                            !! Runtime parameter `DRAG_BG_VEL`.
  real    :: BBL_thick_min  !< The minimum bottom boundary layer thickness [H ~> m or kg m-2].
                            !! This might be Kv / (cdrag * drag_bg_vel) to give
                            !! Kv as the minimum near-bottom viscosity.
  real    :: Htbl_shelf     !< A nominal thickness of the surface boundary layer for use
                            !! in calculating the near-surface velocity [H ~> m or kg m-2].
  real    :: Htbl_shelf_min !< The minimum surface boundary layer thickness [H ~> m or kg m-2].
  real    :: KV_BBL_min     !< The minimum viscosity in the bottom boundary layer [Z2 T-1 ~> m2 s-1].
  real    :: KV_TBL_min     !< The minimum viscosity in the top boundary layer [Z2 T-1 ~> m2 s-1].
  logical :: bottomdraglaw  !< If true, the  bottom stress is calculated with a
                            !! drag law c_drag*|u|*u. The velocity magnitude
                            !! may be an assumed value or it may be based on the
                            !! actual velocity in the bottommost `HBBL`, depending
                            !! on whether linear_drag is true.
                            !! Runtime parameter `BOTTOMDRAGLAW`.
  logical :: BBL_use_EOS    !< If true, use the equation of state in determining
                            !! the properties of the bottom boundary layer.
  logical :: linear_drag    !< If true, the drag law is cdrag*`DRAG_BG_VEL`*u.
                            !! Runtime parameter `LINEAR_DRAG`.
  logical :: Channel_drag   !< If true, the drag is exerted directly on each
                            !! layer according to what fraction of the bottom
                            !! they overlie.
  logical :: correct_BBL_bounds !< If true, uses the correct bounds on the BBL thickness and
                            !! viscosity so that the bottom layer feels the intended drag.
  logical :: RiNo_mix       !< If true, use Richardson number dependent mixing.
  logical :: dynamic_viscous_ML !< If true, use a bulk Richardson number criterion to
                            !! determine the mixed layer thickness for viscosity.
  real    :: bulk_Ri_ML     !< The bulk mixed layer used to determine the
                            !! thickness of the viscous mixed layer.  Nondim.
  real    :: omega          !<   The Earth's rotation rate [T-1 ~> s-1].
  real    :: ustar_min      !< A minimum value of ustar to avoid numerical
                            !! problems [Z T-1 ~> m s-1].  If the value is small enough,
                            !! this should not affect the solution.
  real    :: TKE_decay      !< The ratio of the natural Ekman depth to the TKE
                            !! decay scale, nondimensional.
  real    :: omega_frac     !<   When setting the decay scale for turbulence, use
                            !! this fraction of the absolute rotation rate blended
                            !! with the local value of f, as sqrt((1-of)*f^2 + of*4*omega^2).
  logical :: answers_2018   !< If true, use the order of arithmetic and expressions that recover the
                            !! answers from the end of 2018.  Otherwise, use updated and more robust
                            !! forms of the same expressions.
  logical :: debug          !< If true, write verbose checksums for debugging purposes.
  logical :: BBL_use_tidal_bg !< If true, use a tidal background amplitude for the bottom velocity
                            !! when computing the bottom stress.
  character(len=200) :: inputdir !< The directory for input files.
  type(ocean_OBC_type), pointer :: OBC => NULL() !< Open boundaries control structure
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                            !! regulate the timing of diagnostic output.
  ! Allocatable data arrays
  real, allocatable, dimension(:,:) :: tideamp !< RMS tidal amplitude at h points [Z T-1 ~> m s-1]
  ! Diagnostic arrays
  real, allocatable, dimension(:,:) :: bbl_u !< BBL mean U current [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:) :: bbl_v !< BBL mean V current [L T-1 ~> m s-1]
  !>@{ Diagnostics handles
  integer :: id_bbl_thick_u = -1, id_kv_bbl_u = -1, id_bbl_u = -1
  integer :: id_bbl_thick_v = -1, id_kv_bbl_v = -1, id_bbl_v = -1
  integer :: id_Ray_u = -1, id_Ray_v = -1
  integer :: id_nkml_visc_u = -1, id_nkml_visc_v = -1
  !>@}
end type set_visc_CS

contains

!> Calculates the thickness of the bottom boundary layer and the viscosity within that layer.
!!
!! A drag law is used, either linearized about an assumed bottom velocity or using the
!! actual near-bottom velocities combined with an assumed unresolved velocity.  The bottom
!! boundary layer thickness is limited by a combination of stratification and rotation, as
!! in the paper of Killworth and Edwards, JPO 1999. It is not necessary to calculate the
!! thickness and viscosity every time step; instead previous values may be used.
!!
!! \section set_viscous_BBL Viscous Bottom Boundary Layer
!!
!! If set_visc_cs.bottomdraglaw is True then a bottom boundary layer viscosity and thickness
!! are calculated so that the bottom stress is
!! \f[
!! \mathbf{\tau}_b = C_d | U_{bbl} | \mathbf{u}_{bbl}
!! \f]
!! If set_visc_cs.bottomdraglaw is True then the term \f$|U_{bbl}|\f$ is set equal to the
!! value in set_visc_cs.drag_bg_vel so that \f$C_d |U_{bbl}|\f$ becomes a Rayleigh bottom drag.
!! Otherwise \f$|U_{bbl}|\f$ is found by averaging the flow over the bottom set_visc_cs.hbbl
!! of the model, adding the amplitude of tides set_visc_cs.tideamp and a constant
!! set_visc_cs.drag_bg_vel. For these calculations the vertical grid at the velocity
!! component locations is found by
!! \f[
!! \begin{array}{ll}
!! \frac{2 h^- h^+}{h^- + h^+} & u \left( h^+ - h^-\right) >= 0
!! \\
!! \frac{1}{2} \left( h^- + h^+ \right) &  u \left( h^+ - h^-\right) < 0
!! \end{array}
!! \f]
!! which biases towards the thin cell if the thin cell is upwind. Biasing the grid toward
!! thin upwind cells helps increase the effect of viscosity and inhibits flow out of these
!! thin cells.
!!
!! After diagnosing \f$|U_{bbl}|\f$ over a fixed depth an active viscous boundary layer
!! thickness is found using the ideas of Killworth and Edwards, 1999 (hereafter KW99).
!! KW99 solve the equation
!! \f[
!! \left( \frac{h_{bbl}}{h_f} \right)^2 + \frac{h_{bbl}}{h_N} = 1
!! \f]
!! for the boundary layer depth \f$h_{bbl}\f$. Here
!! \f[
!! h_f = \frac{C_n u_*}{f}
!! \f]
!! is the rotation controlled boundary layer depth in the absence of stratification.
!! \f$u_*\f$ is the surface friction speed given by
!! \f[
!! u_*^2 = C_d |U_{bbl}|^2
!! \f]
!! and is a function of near bottom model flow.
!! \f[
!! h_N = \frac{C_i u_*}{N} = \frac{ (C_i u_* )^2 }{g^\prime}
!! \f]
!! is the stratification controlled boundary layer depth. The non-dimensional parameters
!! \f$C_n=0.5\f$ and \f$C_i=20\f$ are suggested by Zilitinkevich and Mironov, 1996.
!!
!! If a Richardson number dependent mixing scheme is being used, as indicated by
!! set_visc_cs.rino_mix, then the boundary layer thickness is bounded to be no larger
!! than a half of set_visc_cs.hbbl .
!!
!! \todo Channel drag needs to be explained
!!
!! A BBL viscosity is calculated so that the no-slip boundary condition in the vertical
!! viscosity solver implies the stress \f$\mathbf{\tau}_b\f$.
!!
!! \subsection set_viscous_BBL_ref References
!!
!! \arg Killworth, P. D., and N. R. Edwards, 1999:
!! A Turbulent Bottom Boundary Layer Code for Use in Numerical Ocean Models.
!! J. Phys. Oceanogr., 29, 1221-1238,
!! <a href="https://doi.org/10.1175/1520-0485(1999)029<1221:ATBBLC>2.0.CO;2"
!! >doi:10.1175/1520-0485(1999)029<1221:ATBBLC>2.0.CO;2</a>
!! \arg Zilitinkevich, S., Mironov, D.V., 1996:
!! A multi-limit formulation for the equilibrium depth of a stably stratified boundary layer.
!! Boundary-Layer Meteorology 81, 325-351.
!! <a href="https://doi.org/10.1007/BF02430334">doi:10.1007/BF02430334</a>
!!
subroutine set_viscous_BBL(u, v, h, tv, visc, G, GV, US, CS, symmetrize)
  type(ocean_grid_type),    intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),    intent(in)    :: tv   !< A structure containing pointers to any
                                                  !! available thermodynamic fields. Absent fields
                                                  !! have NULL ptrs..
  type(vertvisc_type),      intent(inout) :: visc !< A structure containing vertical viscosities and
                                                  !! related fields.
  type(set_visc_CS),        pointer       :: CS   !< The control structure returned by a previous
                                                  !! call to set_visc_init.
  logical,        optional, intent(in)    :: symmetrize !< If present and true, do extra calculations
                                                  !! of those values in visc that would be
                                                  !! calculated with symmetric memory.

  ! Local variables
  real, dimension(SZIB_(G)) :: &
    ustar, &    !   The bottom friction velocity [Z T-1 ~> m s-1].
    T_EOS, &    !   The temperature used to calculate the partial derivatives
                ! of density with T and S [degC].
    S_EOS, &    !   The salinity used to calculate the partial derivatives
                ! of density with T and S [ppt].
    dR_dT, &    !   Partial derivative of the density in the bottom boundary
                ! layer with temperature [R degC-1 ~> kg m-3 degC-1].
    dR_dS, &    !   Partial derivative of the density in the bottom boundary
                ! layer with salinity [R ppt-1 ~> kg m-3 ppt-1].
    press       !   The pressure at which dR_dT and dR_dS are evaluated [R L2 T-2 ~> Pa].
  real :: htot      ! Sum of the layer thicknesses up to some point [H ~> m or kg m-2].
  real :: htot_vel  ! Sum of the layer thicknesses up to some point [H ~> m or kg m-2].

  real :: Rhtot ! Running sum of thicknesses times the layer potential
                ! densities [H R ~> kg m-2 or kg2 m-5].
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    D_u, &      ! Bottom depth interpolated to u points [Z ~> m].
    mask_u      ! A mask that disables any contributions from u points that
                ! are land or past open boundary conditions [nondim], 0 or 1.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    D_v, &      ! Bottom depth interpolated to v points [Z ~> m].
    mask_v      ! A mask that disables any contributions from v points that
                ! are land or past open boundary conditions [nondim], 0 or 1.
  real, dimension(SZIB_(G),SZK_(GV)) :: &
    h_at_vel, & ! Layer thickness at a velocity point, using an upwind-biased
                ! second order accurate estimate based on the previous velocity
                ! direction [H ~> m or kg m-2].
    h_vel, &    ! Arithmetic mean of the layer thicknesses adjacent to a
                ! velocity point [H ~> m or kg m-2].
    T_vel, &    ! Arithmetic mean of the layer temperatures adjacent to a
                ! velocity point [degC].
    S_vel, &    ! Arithmetic mean of the layer salinities adjacent to a
                ! velocity point [ppt].
    Rml_vel     ! Arithmetic mean of the layer coordinate densities adjacent
                ! to a velocity point [R ~> kg m-3].

  real :: h_vel_pos        ! The arithmetic mean thickness at a velocity point
                           ! plus H_neglect to avoid 0 values [H ~> m or kg m-2].
  real :: ustarsq          ! 400 times the square of ustar, times
                           ! Rho0 divided by G_Earth and the conversion
                           ! from m to thickness units [H R ~> kg m-2 or kg2 m-5].
  real :: cdrag_sqrt_Z     ! Square root of the drag coefficient, times a unit conversion
                           ! factor from lateral lengths to vertical depths [Z L-1 ~> 1].
  real :: cdrag_sqrt       ! Square root of the drag coefficient [nondim].
  real :: oldfn            ! The integrated energy required to
                           ! entrain up to the bottom of the layer,
                           ! divided by G_Earth [H R ~> kg m-2 or kg2 m-5].
  real :: Dfn              ! The increment in oldfn for entraining
                           ! the layer [H R ~> kg m-2 or kg2 m-5].
  real :: Dh               ! The increment in layer thickness from
                           ! the present layer [H ~> m or kg m-2].
  real :: bbl_thick        ! The thickness of the bottom boundary layer [H ~> m or kg m-2].
  real :: bbl_thick_Z      ! The thickness of the bottom boundary layer [Z ~> m].
  real :: kv_bbl           ! The bottom boundary layer viscosity [Z2 T-1 ~> m2 s-1].
  real :: C2f              ! C2f = 2*f at velocity points [T-1 ~> s-1].

  real :: U_bg_sq          ! The square of an assumed background
                           ! velocity, for calculating the mean
                           ! magnitude near the bottom for use in the
                           ! quadratic bottom drag [L2 T-2 ~> m2 s-2].
  real :: hwtot            ! Sum of the thicknesses used to calculate
                           ! the near-bottom velocity magnitude [H ~> m or kg m-2].
  real :: hutot            ! Running sum of thicknesses times the
                           ! velocity magnitudes [H T T-1 ~> m2 s-1 or kg m-1 s-1].
  real :: Thtot            ! Running sum of thickness times temperature [degC H ~> degC m or degC kg m-2].
  real :: Shtot            ! Running sum of thickness times salinity [ppt H ~> ppt m or ppt kg m-2].
  real :: hweight          ! The thickness of a layer that is within Hbbl
                           ! of the bottom [H ~> m or kg m-2].
  real :: v_at_u, u_at_v   ! v at a u point or vice versa [L T-1 ~> m s-1].
  real :: Rho0x400_G       ! 400*Rho0/G_Earth, times unit conversion factors
                           ! [R T2 H Z-2 ~> kg s2 m-4 or kg2 s2 m-7].
                           ! The 400 is a constant proposed by Killworth and Edwards, 1999.
  real, dimension(SZI_(G),SZJ_(G),max(GV%nk_rho_varies,1)) :: &
    Rml                    ! The mixed layer coordinate density [R ~> kg m-3].
  real :: p_ref(SZI_(G))   !   The pressure used to calculate the coordinate
                           ! density [R L2 T-2 ~> Pa] (usually set to 2e7 Pa = 2000 dbar).

  real :: D_vel            ! The bottom depth at a velocity point [H ~> m or kg m-2].
  real :: Dp, Dm           ! The depths at the edges of a velocity cell [H ~> m or kg m-2].
  real :: a                ! a is the curvature of the bottom depth across a
                           ! cell, times the cell width squared [H ~> m or kg m-2].
  real :: a_3, a_12        ! a/3 and a/12 [H ~> m or kg m-2].
  real :: C24_a            ! 24/a [H-1 ~> m-1 or m2 kg-1].
  real :: slope            ! The absolute value of the bottom depth slope across
                           ! a cell times the cell width [H ~> m or kg m-2].
  real :: apb_4a, ax2_3apb ! Various nondimensional ratios of a and slope.
  real :: a2x48_apb3, Iapb, Ibma_2 ! Combinations of a and slope [H-1 ~> m-1 or m2 kg-1].
  ! All of the following "volumes" have units of thickness because they are normalized
  ! by the full horizontal area of a velocity cell.
  real :: Vol_open         ! The cell volume above which it is open [H ~> m or kg m-2].
  real :: Vol_direct       ! With less than Vol_direct [H ~> m or kg m-2], there is a direct
                           ! solution of a cubic equation for L.
  real :: Vol_2_reg        ! The cell volume above which there are two separate
                           ! open areas that must be integrated [H ~> m or kg m-2].
  real :: vol              ! The volume below the interface whose normalized
                           ! width is being sought [H ~> m or kg m-2].
  real :: vol_below        ! The volume below the interface below the one that
                           ! is currently under consideration [H ~> m or kg m-2].
  real :: Vol_err          ! The error in the volume with the latest estimate of
                           ! L, or the error for the interface below [H ~> m or kg m-2].
  real :: Vol_quit         ! The volume error below which to quit iterating [H ~> m or kg m-2].
  real :: Vol_tol          ! A volume error tolerance [H ~> m or kg m-2].
  real :: L(SZK_(GV)+1)    ! The fraction of the full cell width that is open at
                           ! the depth of each interface [nondim].
  real :: L_direct         ! The value of L above volume Vol_direct [nondim].
  real :: L_max, L_min     ! Upper and lower bounds on the correct value for L  [nondim].
  real :: Vol_err_max      ! The volume errors for the upper and lower bounds on
  real :: Vol_err_min      ! the correct value for L [H ~> m or kg m-2].
  real :: Vol_0            ! A deeper volume with known width L0 [H ~> m or kg m-2].
  real :: L0               ! The value of L above volume Vol_0 [nondim].
  real :: dVol             ! vol - Vol_0 [H ~> m or kg m-2].
  real :: dV_dL2           ! The partial derivative of volume with L squared
                           ! evaluated at L=L0 [H ~> m or kg m-2].
  real :: h_neglect        ! A thickness that is so small it is usually lost
                           ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: ustH             ! ustar converted to units of H T-1 [H T-1 ~> m s-1 or kg m-2 s-1].
  real :: root             ! A temporary variable [H T-1 ~> m s-1 or kg m-2 s-1].

  real :: Cell_width       ! The transverse width of the velocity cell [L ~> m].
  real :: Rayleigh         ! A nondimensional value that is multiplied by the layer's
                           ! velocity magnitude to give the Rayleigh drag velocity, times
                           ! a lateral to vertical distance conversion factor [Z L-1 ~> 1].
  real :: gam              ! The ratio of the change in the open interface width
                           ! to the open interface width atop a cell [nondim].
  real :: BBL_frac         ! The fraction of a layer's drag that goes into the
                           ! viscous bottom boundary layer [nondim].
  real :: BBL_visc_frac    ! The fraction of all the drag that is expressed as
                           ! a viscous bottom boundary layer [nondim].
  real, parameter :: C1_3 = 1.0/3.0, C1_6 = 1.0/6.0, C1_12 = 1.0/12.0
  real :: C2pi_3           ! An irrational constant, 2/3 pi.
  real :: tmp              ! A temporary variable.
  real :: tmp_val_m1_to_p1
  real :: curv_tol         ! Numerator of curvature cubed, used to estimate
                           ! accuracy of a single L(:) Newton iteration
  logical :: use_L0, do_one_L_iter    ! Control flags for L(:) Newton iteration
  logical :: use_BBL_EOS, do_i(SZIB_(G))
  integer, dimension(2) :: EOSdom ! The computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, m, n, K2, nkmb, nkml
  integer :: itt, maxitt=20
  type(ocean_OBC_type), pointer :: OBC => NULL()

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nkmb = GV%nk_rho_varies ; nkml = GV%nkml
  h_neglect = GV%H_subroundoff
  Rho0x400_G = 400.0*(GV%Rho0 / (US%L_to_Z**2 * GV%g_Earth)) * GV%Z_to_H
  Vol_quit = 0.9*GV%Angstrom_H + h_neglect
  C2pi_3 = 8.0*atan(1.0)/3.0

  if (.not.associated(CS)) call MOM_error(FATAL,"MOM_set_viscosity(BBL): "//&
         "Module must be initialized before it is used.")
  if (.not.CS%bottomdraglaw) return

  if (present(symmetrize)) then ; if (symmetrize) then
    Jsq = js-1 ; Isq = is-1
  endif ; endif

  if (CS%debug) then
    call uvchksum("Start set_viscous_BBL [uv]", u, v, G%HI, haloshift=1, scale=US%L_T_to_m_s)
    call hchksum(h,"Start set_viscous_BBL h", G%HI, haloshift=1, scale=GV%H_to_m)
    if (associated(tv%T)) call hchksum(tv%T, "Start set_viscous_BBL T", G%HI, haloshift=1)
    if (associated(tv%S)) call hchksum(tv%S, "Start set_viscous_BBL S", G%HI, haloshift=1)
  endif

  use_BBL_EOS = associated(tv%eqn_of_state) .and. CS%BBL_use_EOS
  OBC => CS%OBC

  U_bg_sq = CS%drag_bg_vel * CS%drag_bg_vel
  cdrag_sqrt = sqrt(CS%cdrag)
  cdrag_sqrt_Z = US%L_to_Z * sqrt(CS%cdrag)
  K2 = max(nkmb+1, 2)

!  With a linear drag law, the friction velocity is already known.
!  if (CS%linear_drag) ustar(:) = cdrag_sqrt_Z*CS%drag_bg_vel

  if ((nkml>0) .and. .not.use_BBL_EOS) then
    EOSdom(1) = Isq - (G%isd-1) ;  EOSdom(2) = G%iec+1 - (G%isd-1)
    do i=Isq,Ieq+1 ; p_ref(i) = tv%P_Ref ; enddo
    !$OMP parallel do default(shared)
    do k=1,nkmb ; do j=Jsq,Jeq+1
      call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_ref, Rml(:,j,k), tv%eqn_of_state, &
                             EOSdom)
    enddo ; enddo
  endif

  !$OMP parallel do default(shared)
  do J=js-1,je ; do i=is-1,ie+1
    D_v(i,J) = 0.5*(G%bathyT(i,j) + G%bathyT(i,j+1))
    mask_v(i,J) = G%mask2dCv(i,J)
  enddo ; enddo
  !$OMP parallel do default(shared)
  do j=js-1,je+1 ; do I=is-1,ie
    D_u(I,j) = 0.5*(G%bathyT(i,j) + G%bathyT(i+1,j))
    mask_u(I,j) = G%mask2dCu(I,j)
  enddo ; enddo

  if (associated(OBC)) then ; do n=1,OBC%number_of_segments
    if (.not. OBC%segment(n)%on_pe) cycle
    ! Use a one-sided projection of bottom depths at OBC points.
    I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
    if (OBC%segment(n)%is_N_or_S .and. (J >= js-1) .and. (J <= je)) then
      do i = max(is-1,OBC%segment(n)%HI%isd), min(ie+1,OBC%segment(n)%HI%ied)
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) D_v(i,J) = G%bathyT(i,j)
        if (OBC%segment(n)%direction == OBC_DIRECTION_S) D_v(i,J) = G%bathyT(i,j+1)
      enddo
    elseif (OBC%segment(n)%is_E_or_W .and. (I >= is-1) .and. (I <= ie)) then
      do j = max(js-1,OBC%segment(n)%HI%jsd), min(je+1,OBC%segment(n)%HI%jed)
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) D_u(I,j) = G%bathyT(i,j)
        if (OBC%segment(n)%direction == OBC_DIRECTION_W) D_u(I,j) = G%bathyT(i+1,j)
      enddo
    endif
  enddo ; endif
  if (associated(OBC)) then ; do n=1,OBC%number_of_segments
    ! Now project bottom depths across cell-corner points in the OBCs.  The two
    ! projections have to occur in sequence and can not be combined easily.
    if (.not. OBC%segment(n)%on_pe) cycle
    ! Use a one-sided projection of bottom depths at OBC points.
    I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
    if (OBC%segment(n)%is_N_or_S .and. (J >= js-1) .and. (J <= je)) then
      do I = max(is-1,OBC%segment(n)%HI%IsdB), min(ie,OBC%segment(n)%HI%IedB)
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          D_u(I,j+1) = D_u(I,j) ;  mask_u(I,j+1) = 0.0
        elseif (OBC%segment(n)%direction == OBC_DIRECTION_S) then
          D_u(I,j) = D_u(I,j+1) ; mask_u(I,j) = 0.0
        endif
      enddo
    elseif (OBC%segment(n)%is_E_or_W .and. (I >= is-1) .and. (I <= ie)) then
      do J = max(js-1,OBC%segment(n)%HI%JsdB), min(je,OBC%segment(n)%HI%JedB)
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          D_v(i+1,J) = D_v(i,J) ; mask_v(i+1,J) = 0.0
        elseif (OBC%segment(n)%direction == OBC_DIRECTION_W) then
          D_v(i,J) = D_v(i+1,J) ;  mask_v(i,J) = 0.0
        endif
      enddo
    endif
  enddo ; endif

  if (.not.use_BBL_EOS) Rml_vel(:,:) = 0.0

  !$OMP parallel do default(private) shared(u,v,h,tv,visc,G,GV,US,CS,Rml,nz,nkmb, &
  !$OMP                                     nkml,Isq,Ieq,Jsq,Jeq,h_neglect,Rho0x400_G,C2pi_3, &
  !$OMP                                     U_bg_sq,cdrag_sqrt_Z,cdrag_sqrt,K2,use_BBL_EOS,   &
  !$OMP                                     OBC,maxitt,D_u,D_v,mask_u,mask_v) &
  !$OMP                              firstprivate(Vol_quit)
  do j=Jsq,Jeq ; do m=1,2

    if (m==1) then
      ! m=1 refers to u-points
      if (j<G%Jsc) cycle
      is = Isq ; ie = Ieq
      do i=is,ie
        do_i(i) = .false.
        if (G%mask2dCu(I,j) > 0) do_i(i) = .true.
      enddo
    else
      ! m=2 refers to v-points
      is = G%isc ; ie = G%iec
      do i=is,ie
        do_i(i) = .false.
        if (G%mask2dCv(i,J) > 0) do_i(i) = .true.
      enddo
    endif

    ! Calculate thickness at velocity points (u or v depending on value of m).
    ! Also interpolate the ML density or T/S properties.
    if (m==1) then ! u-points
      do k=1,nz ; do I=is,ie
        if (do_i(I)) then
          if (u(I,j,k) *(h(i+1,j,k) - h(i,j,k)) >= 0) then
            ! If the flow is from thin to thick then bias towards the thinner thickness
            h_at_vel(I,k) = 2.0*h(i,j,k)*h(i+1,j,k) / &
                            (h(i,j,k) + h(i+1,j,k) + h_neglect)
          else
            ! If the flow is from thick to thin then use the simple average thickness
            h_at_vel(I,k) = 0.5 * (h(i,j,k) + h(i+1,j,k))
          endif
        endif
        h_vel(I,k) = 0.5 * (h(i,j,k) + h(i+1,j,k))
      enddo ; enddo
      if (use_BBL_EOS) then ; do k=1,nz ; do I=is,ie
        ! Perhaps these should be thickness weighted.
        T_vel(I,k) = 0.5 * (tv%T(i,j,k) + tv%T(i+1,j,k))
        S_vel(I,k) = 0.5 * (tv%S(i,j,k) + tv%S(i+1,j,k))
      enddo ; enddo ; else ; do k=1,nkmb ; do I=is,ie
        Rml_vel(I,k) = 0.5 * (Rml(i,j,k) + Rml(i+1,j,k))
      enddo ; enddo ; endif
    else ! v-points
      do k=1,nz ; do i=is,ie
        if (do_i(i)) then
          if (v(i,J,k) * (h(i,j+1,k) - h(i,j,k)) >= 0) then
            ! If the flow is from thin to thick then bias towards the thinner thickness
            h_at_vel(i,k) = 2.0*h(i,j,k)*h(i,j+1,k) / &
                            (h(i,j,k) + h(i,j+1,k) + h_neglect)
          else
            ! If the flow is from thick to thin then use the simple average thickness
            h_at_vel(i,k) = 0.5 * (h(i,j,k) + h(i,j+1,k))
          endif
        endif
        h_vel(i,k) = 0.5 * (h(i,j,k) + h(i,j+1,k))
      enddo ; enddo
      if (use_BBL_EOS) then ; do k=1,nz ; do i=is,ie
        T_vel(i,k) = 0.5 * (tv%T(i,j,k) + tv%T(i,j+1,k))
        S_vel(i,k) = 0.5 * (tv%S(i,j,k) + tv%S(i,j+1,k))
      enddo ; enddo ; else ; do k=1,nkmb ; do i=is,ie
        Rml_vel(i,k) = 0.5 * (Rml(i,j,k) + Rml(i,j+1,k))
      enddo ; enddo ; endif
    endif

    if (associated(OBC)) then ; if (OBC%number_of_segments > 0) then
      ! Apply a zero gradient projection of thickness across OBC points.
      if (m==1) then
        do I=is,ie ; if (do_i(I) .and. (OBC%segnum_u(I,j) /= OBC_NONE)) then
          if (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_E) then
            do k=1,nz
              h_at_vel(I,k) = h(i,j,k) ; h_vel(I,k) = h(i,j,k)
            enddo
            if (use_BBL_EOS) then
              do k=1,nz
                T_vel(I,k) = tv%T(i,j,k) ; S_vel(I,k) = tv%S(i,j,k)
              enddo
            else
              do k=1,nkmb
                Rml_vel(I,k) = Rml(i,j,k)
              enddo
            endif
          elseif (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_W) then
            do k=1,nz
              h_at_vel(I,k) = h(i+1,j,k) ; h_vel(I,k) = h(i+1,j,k)
            enddo
            if (use_BBL_EOS) then
              do k=1,nz
                T_vel(I,k) = tv%T(i+1,j,k) ; S_vel(I,k) = tv%S(i+1,j,k)
              enddo
            else
              do k=1,nkmb
                Rml_vel(I,k) = Rml(i+1,j,k)
              enddo
            endif
          endif
        endif ; enddo
      else
        do i=is,ie ; if (do_i(i) .and. (OBC%segnum_v(i,J) /= OBC_NONE)) then
          if (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_N) then
            do k=1,nz
              h_at_vel(i,k) = h(i,j,k) ; h_vel(i,k) = h(i,j,k)
            enddo
            if (use_BBL_EOS) then
              do k=1,nz
                T_vel(i,k) = tv%T(i,j,k) ; S_vel(i,k) = tv%S(i,j,k)
              enddo
            else
              do k=1,nkmb
                Rml_vel(i,k) = Rml(i,j,k)
              enddo
            endif
          elseif (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_S) then
            do k=1,nz
              h_at_vel(i,k) = h(i,j+1,k) ; h_vel(i,k) = h(i,j+1,k)
            enddo
            if (use_BBL_EOS) then
              do k=1,nz
                T_vel(i,k) = tv%T(i,j+1,k) ; S_vel(i,k) = tv%S(i,j+1,k)
              enddo
            else
              do k=1,nkmb
                Rml_vel(i,k) = Rml(i,j+1,k)
              enddo
            endif
          endif
        endif ; enddo
      endif
    endif ; endif

    if (use_BBL_EOS .or. .not.CS%linear_drag) then
      ! Calculate the mean velocity magnitude over the bottommost CS%Hbbl of
      ! the water column for determining the quadratic bottom drag.
      ! Used in ustar(i)
      do i=is,ie ; if (do_i(i)) then
        htot_vel = 0.0 ; hwtot = 0.0 ; hutot = 0.0
        Thtot = 0.0 ; Shtot = 0.0
        do k=nz,1,-1

          if (htot_vel>=CS%Hbbl) exit ! terminate the k loop

          hweight = MIN(CS%Hbbl - htot_vel, h_at_vel(i,k))
          if (hweight < 1.5*GV%Angstrom_H + h_neglect) cycle

          htot_vel  = htot_vel + h_at_vel(i,k)
          hwtot = hwtot + hweight

          if ((.not.CS%linear_drag) .and. (hweight >= 0.0)) then ; if (m==1) then
            v_at_u = set_v_at_u(v, h, G, GV, i, j, k, mask_v, OBC)
            if (CS%BBL_use_tidal_bg) then
              U_bg_sq = 0.5*( G%mask2dT(i,j)*(CS%tideamp(i,j)*CS%tideamp(i,j))+ &
                              G%mask2dT(i+1,j)*(CS%tideamp(i+1,j)*CS%tideamp(i+1,j)) )
            endif
            hutot = hutot + hweight * sqrt(u(I,j,k)*u(I,j,k) + &
                                           v_at_u*v_at_u + U_bg_sq)
          else
            u_at_v = set_u_at_v(u, h, G, GV, i, j, k, mask_u, OBC)
            if (CS%BBL_use_tidal_bg) then
              U_bg_sq = 0.5*( G%mask2dT(i,j)*(CS%tideamp(i,j)*CS%tideamp(i,j))+ &
                              G%mask2dT(i,j+1)*(CS%tideamp(i,j+1)*CS%tideamp(i,j+1)) )
            endif
            hutot = hutot + hweight * sqrt(v(i,J,k)*v(i,J,k) + &
                                           u_at_v*u_at_v + U_bg_sq)
          endif ; endif

          if (use_BBL_EOS .and. (hweight >= 0.0)) then
            Thtot = Thtot + hweight * T_vel(i,k)
            Shtot = Shtot + hweight * S_vel(i,k)
          endif
        enddo ! end of k loop

        ! Set u* based on u*^2 = Cdrag u_bbl^2
        if (.not.CS%linear_drag .and. (hwtot > 0.0)) then
          ustar(i) = cdrag_sqrt_Z*hutot/hwtot
        else
          ustar(i) = cdrag_sqrt_Z*CS%drag_bg_vel
        endif

        if (use_BBL_EOS) then ; if (hwtot > 0.0) then
          T_EOS(i) = Thtot/hwtot ; S_EOS(i) = Shtot/hwtot
        else
          T_EOS(i) = 0.0 ; S_EOS(i) = 0.0
        endif ; endif

        ! Diagnostic BBL flow speed at u- and v-points.
        if (CS%id_bbl_u>0 .and. m==1) then
          if (hwtot > 0.0) CS%bbl_u(I,j) = hutot/hwtot
        elseif (CS%id_bbl_v>0 .and. m==2) then
          if (hwtot > 0.0) CS%bbl_v(i,J) = hutot/hwtot
        endif

      endif ; enddo
    else
      do i=is,ie ; ustar(i) = cdrag_sqrt_Z*CS%drag_bg_vel ; enddo
    endif ! Not linear_drag

    if (use_BBL_EOS) then
      if (associated(tv%p_surf)) then
        if (m==1) then ; do i=is,ie ; press(I) = 0.5*(tv%p_surf(i,j) + tv%p_surf(i+1,j)) ; enddo
        else ; do i=is,ie ; press(i) = 0.5*(tv%p_surf(i,j) + tv%p_surf(i,j+1)) ; enddo ; endif
      else
        do i=is,ie ; press(i) = 0.0 ; enddo
      endif
      do i=is,ie ; if (.not.do_i(i)) then ; T_EOS(i) = 0.0 ; S_EOS(i) = 0.0 ; endif ; enddo
      do k=1,nz ; do i=is,ie
        press(i) = press(i) + (GV%H_to_RZ*GV%g_Earth) * h_vel(i,k)
      enddo ; enddo
      call calculate_density_derivs(T_EOS, S_EOS, press, dR_dT, dR_dS, tv%eqn_of_state, &
                                    (/is-G%IsdB+1,ie-G%IsdB+1/) )
    endif

    ! Find a BBL thickness given by equation 2.20 of Killworth and Edwards, 1999:
    !    ( f h / Cn u* )^2 + ( N h / Ci u* ) = 1
    ! where Cn=0.5 and Ci=20 (constants suggested by Zilitinkevich and Mironov, 1996).
    ! Eq. 2.20 can be expressed in terms of boundary layer thicknesses limited by
    ! rotation (h_f) and stratification (h_N):
    !    ( h / h_f )^2 + ( h / h_N ) = 1
    ! When stratification dominates h_N<<h_f, and vice versa.
    do i=is,ie ; if (do_i(i)) then
      ! The 400.0 in this expression is the square of a Ci introduced in KW99, eq. 2.22.
      ustarsq = Rho0x400_G * ustar(i)**2 ! Note not in units of u*^2 but [H R]
      htot = 0.0

      ! Calculate the thickness of a stratification limited BBL ignoring rotation:
      !   h_N = Ci u* / N          (limit of KW99 eq. 2.20 for |f|->0)
      ! For layer mode, N^2 = g'/h. Since (Ci u*)^2 = (h_N N)^2 = h_N g' then
      !   h_N = (Ci u*)^2 / g'     (KW99, eq, 2.22)
      ! Starting from the bottom, integrate the stratification upward until h_N N balances Ci u*
      ! or in layer mode
      !   h_N Delta rho ~ (Ci u*)^2 rho0 / g
      ! where the rhs is stored in variable ustarsq.
      ! The method was described in Stephens and Hallberg 2000 (unpublished and lost manuscript).
      if (use_BBL_EOS) then
        Thtot = 0.0 ; Shtot = 0.0 ; oldfn = 0.0
        do k=nz,2,-1
          if (h_at_vel(i,k) <= 0.0) cycle

          ! Delta rho * h_bbl assuming everything below is homogenized
          oldfn = dR_dT(i)*(Thtot - T_vel(i,k)*htot) + &
                  dR_dS(i)*(Shtot - S_vel(i,k)*htot)
          if (oldfn >= ustarsq) exit

          ! Local Delta rho * h_bbl  at interface
          Dfn = (dR_dT(i)*(T_vel(i,k) - T_vel(i,k-1)) + &
                 dR_dS(i)*(S_vel(i,k) - S_vel(i,k-1))) * &
                (h_at_vel(i,k) + htot)

          if ((oldfn + Dfn) <= ustarsq) then
            ! Use whole layer
            Dh = h_at_vel(i,k)
          else
            ! Use only part of the layer
            Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn) / (Dfn))
          endif

          ! Increment total BBL thickness and cumulative T and S
          htot = htot + Dh
          Thtot = Thtot + T_vel(i,k)*Dh ; Shtot = Shtot + S_vel(i,k)*Dh
        enddo
        if ((oldfn < ustarsq) .and. h_at_vel(i,1) > 0.0) then
          ! Layer 1 might be part of the BBL.
          if (dR_dT(i) * (Thtot - T_vel(i,1)*htot) + &
              dR_dS(i) * (Shtot - S_vel(i,1)*htot) < ustarsq) &
            htot = htot + h_at_vel(i,1)
        endif ! Examination of layer 1.
      else  ! Use Rlay and/or the coordinate density as density variables.
        Rhtot = 0.0
        do k=nz,K2,-1
          oldfn = Rhtot - GV%Rlay(k)*htot
          Dfn = (GV%Rlay(k) - GV%Rlay(k-1))*(h_at_vel(i,k)+htot)

          if (oldfn >= ustarsq) then
            cycle
          elseif ((oldfn + Dfn) <= ustarsq) then
            Dh = h_at_vel(i,k)
          else
            Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn) / (Dfn))
          endif

          htot = htot + Dh
          Rhtot = Rhtot + GV%Rlay(k)*Dh
        enddo
        if (nkml>0) then
          do k=nkmb,2,-1
            oldfn = Rhtot - Rml_vel(i,k)*htot
            Dfn = (Rml_vel(i,k) - Rml_vel(i,k-1)) * (h_at_vel(i,k)+htot)

            if (oldfn >= ustarsq) then
              cycle
            elseif ((oldfn + Dfn) <= ustarsq) then
              Dh = h_at_vel(i,k)
            else
              Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn) / (Dfn))
            endif

            htot = htot + Dh
            Rhtot = Rhtot + Rml_vel(i,k)*Dh
          enddo
          if (Rhtot - Rml_vel(i,1)*htot < ustarsq) htot = htot + h_at_vel(i,1)
        else
          if (Rhtot - GV%Rlay(1)*htot < ustarsq) htot = htot + h_at_vel(i,1)
        endif
      endif ! use_BBL_EOS

      ! Value of 2*f at u- or v-points.
      if (m==1) then ; C2f = G%CoriolisBu(I,J-1) + G%CoriolisBu(I,J)
      else ; C2f = G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J) ; endif

      ! The thickness of a rotation limited BBL ignoring stratification is
      !   h_f ~ Cn u* / f        (limit of KW99 eq. 2.20 for N->0).
      ! The buoyancy limit of BBL thickness (h_N) is already in the variable htot from above.
      ! Substituting x = h_N/h into KW99 eq. 2.20 yields the quadratic
      !   x^2 - x = (h_N / h_f)^2
      ! for which the positive root is
      !   xp = 1/2 + sqrt( 1/4 + (h_N/h_f)^2 )
      ! and thus h_bbl = h_N / xp . Since h_f = Cn u*/f and Cn=0.5
      !   xp = 1/2 + sqrt( 1/4 + (2 f h_N/u*)^2 )
      ! To avoid dividing by zero if u*=0 then
      !   xp u* = 1/2 u* + sqrt( 1/4 u*^2 + (2 f h_N)^2 )
      if (CS%cdrag * U_bg_sq <= 0.0) then
        ! This avoids NaNs and overflows, and could be used in all cases,
        ! but is not bitwise identical to the current code.
        ustH = ustar(i)*GV%Z_to_H ; root = sqrt(0.25*ustH**2 + (htot*C2f)**2)
        if (htot*ustH <= (CS%BBL_thick_min+h_neglect) * (0.5*ustH + root)) then
          bbl_thick = CS%BBL_thick_min
        else
          ! The following expression reads
          !   h_bbl = h_N u* / ( 1/2 u* + sqrt( 1/4 u*^2 + ( 2 f h_N )^2 ) )
          ! which is h_bbl = h_N u*/(xp u*) as described above.
          bbl_thick = (htot * ustH) / (0.5*ustH + root)
        endif
      else
        ! The following expression reads
        !   h_bbl = h_N / ( 1/2 + sqrt( 1/4 + ( 2 f h_N / u* )^2 ) )
        ! which is h_bbl = h_N/xp as described above.
        bbl_thick = htot / (0.5 + sqrt(0.25 + htot*htot*C2f*C2f/ &
                                       ((ustar(i)*ustar(i)) * (GV%Z_to_H**2)) ) )

        if (bbl_thick < CS%BBL_thick_min) bbl_thick = CS%BBL_thick_min
      endif

      ! If there is Richardson number dependent mixing, that determines
      ! the vertical extent of the bottom boundary layer, and there is no
      ! need to set that scale here.  In fact, viscously reducing the
      ! shears over an excessively large region reduces the efficacy of
      ! the Richardson number dependent mixing.
      ! In other words, if using RiNo_mix then CS%Hbbl acts as an upper bound on
      ! bbl_thick.
      if ((bbl_thick > 0.5*CS%Hbbl) .and. (CS%RiNo_mix)) bbl_thick = 0.5*CS%Hbbl

      if (CS%Channel_drag) then
        ! The drag within the bottommost bbl_thick is applied as a part of
        ! an enhanced bottom viscosity, while above this the drag is applied
        ! directly to the layers in question as a Rayleigh drag term.
        if (m==1) then
          D_vel = D_u(I,j)
          tmp = G%mask2dCu(I,j+1) * D_u(I,j+1)
          Dp = 2.0 * D_vel * tmp / (D_vel + tmp)
          tmp = G%mask2dCu(I,j-1) * D_u(I,j-1)
          Dm = 2.0 * D_vel * tmp / (D_vel + tmp)
        else
          D_vel = D_v(i,J)
          tmp = G%mask2dCv(i+1,J) * D_v(i+1,J)
          Dp = 2.0 * D_vel * tmp / (D_vel + tmp)
          tmp = G%mask2dCv(i-1,J) * D_v(i-1,J)
          Dm = 2.0 * D_vel * tmp / (D_vel + tmp)
        endif
        if (Dm > Dp) then ; tmp = Dp ; Dp = Dm ; Dm = tmp ; endif

        ! Convert the D's to the units of thickness.
        Dp = GV%Z_to_H*Dp ; Dm = GV%Z_to_H*Dm ; D_vel = GV%Z_to_H*D_vel

        a_3 = (Dp + Dm - 2.0*D_vel) ; a = 3.0*a_3 ; a_12 = 0.25*a_3
        slope = Dp - Dm
        ! If the curvature is small enough, there is no reason not to assume
        ! a uniformly sloping or flat bottom.
        if (abs(a) < 1e-2*(slope + CS%BBL_thick_min)) a = 0.0
        ! Each cell extends from x=-1/2 to 1/2, and has a topography
        ! given by D(x) = a*x^2 + b*x + D - a/12.

        ! Calculate the volume above which the entire cell is open and the
        ! other volumes at which the equation that is solved for L changes.
        if (a > 0.0) then
          if (slope >= a) then
            Vol_open = D_vel - Dm ; Vol_2_reg = Vol_open
          else
            tmp = slope/a
            Vol_open = 0.25*slope*tmp + C1_12*a
            Vol_2_reg = 0.5*tmp**2 * (a - C1_3*slope)
          endif
          ! Define some combinations of a & b for later use.
          C24_a = 24.0/a ; Iapb = 1.0/(a+slope)
          apb_4a = (slope+a)/(4.0*a) ; a2x48_apb3 = (48.0*(a*a))*(Iapb**3)
          ax2_3apb = 2.0*C1_3*a*Iapb
        elseif (a == 0.0) then
          Vol_open = 0.5*slope
          if (slope > 0) Iapb = 1.0/slope
        else ! a < 0.0
          Vol_open = D_vel - Dm
          if (slope >= -a) then
            Iapb = 1.0e30 ; if (slope+a /= 0.0) Iapb = 1.0/(a+slope)
            Vol_direct = 0.0 ; L_direct = 0.0 ; C24_a = 0.0
          else
            C24_a = 24.0/a ; Iapb = 1.0/(a+slope)
            L_direct = 1.0 + slope/a ! L_direct < 1 because a < 0
            Vol_direct = -C1_6*a*L_direct**3
          endif
          Ibma_2 = 2.0 / (slope - a)
        endif

        L(nz+1) = 0.0 ; vol = 0.0 ; Vol_err = 0.0 ; BBL_visc_frac = 0.0
        ! Determine the normalized open length at each interface.
        do K=nz,1,-1
          vol_below = vol

          vol = vol + h_vel(i,k)
          h_vel_pos = h_vel(i,k) + h_neglect

          if (vol >= Vol_open) then ; L(K) = 1.0
          elseif (a == 0) then ! The bottom has no curvature.
            L(K) = sqrt(2.0*vol*Iapb)
          elseif (a > 0) then
            ! There may be a minimum depth, and there are
            ! analytic expressions for L for all cases.
            if (vol < Vol_2_reg) then
              ! In this case, there is a contiguous open region and
              !   vol = 0.5*L^2*(slope + a/3*(3-4L)).
              if (a2x48_apb3*vol < 1e-8) then ! Could be 1e-7?
                ! There is a very good approximation here for massless layers.
                L0 = sqrt(2.0*vol*Iapb) ; L(K) = L0*(1.0 + ax2_3apb*L0)
              else
                L(K) = apb_4a * (1.0 - &
                         2.0 * cos(C1_3*acos(a2x48_apb3*vol - 1.0) - C2pi_3))
              endif
              ! To check the answers.
              ! Vol_err = 0.5*(L(K)*L(K))*(slope + a_3*(3.0-4.0*L(K))) - vol
            else ! There are two separate open regions.
              !   vol = slope^2/4a + a/12 - (a/12)*(1-L)^2*(1+2L)
              ! At the deepest volume, L = slope/a, at the top L = 1.
              !L(K) = 0.5 - cos(C1_3*acos(1.0 - C24_a*(Vol_open - vol)) - C2pi_3)
              tmp_val_m1_to_p1 = 1.0 - C24_a*(Vol_open - vol)
              tmp_val_m1_to_p1 = max(-1., min(1., tmp_val_m1_to_p1))
              L(K) = 0.5 - cos(C1_3*acos(tmp_val_m1_to_p1) - C2pi_3)
              ! To check the answers.
              ! Vol_err = Vol_open - a_12*(1.0+2.0*L(K)) * (1.0-L(K))**2 - vol
            endif
          else ! a < 0.
            if (vol <= Vol_direct) then
              ! Both edges of the cell are bounded by walls.
              L(K) = (-0.25*C24_a*vol)**C1_3
            else
              ! x_R is at 1/2 but x_L is in the interior & L is found by solving
              !   vol = 0.5*L^2*(slope + a/3*(3-4L))

              !  Vol_err = 0.5*(L(K+1)*L(K+1))*(slope + a_3*(3.0-4.0*L(K+1))) - vol_below
              ! Change to ...
              !   if (min(Vol_below + Vol_err, vol) <= Vol_direct) then ?
              if (vol_below + Vol_err <= Vol_direct) then
                L0 = L_direct ; Vol_0 = Vol_direct
              else
                L0 = L(K+1) ; Vol_0 = Vol_below + Vol_err
                ! Change to   Vol_0 = min(Vol_below + Vol_err, vol) ?
              endif

              !   Try a relatively simple solution that usually works well
              ! for massless layers.
              dV_dL2 = 0.5*(slope+a) - a*L0 ; dVol = (vol-Vol_0)
           !  dV_dL2 = 0.5*(slope+a) - a*L0 ; dVol = max(vol-Vol_0, 0.0)

              use_L0 = .false.
              do_one_L_iter = .false.
              if (CS%answers_2018) then
                curv_tol = GV%Angstrom_H*dV_dL2**2 &
                           * (0.25 * dV_dL2 * GV%Angstrom_H - a * L0 * dVol)
                do_one_L_iter = (a * a * dVol**3) < curv_tol
              else
                ! The following code is more robust when GV%Angstrom_H=0, but
                ! it changes answers.
                use_L0 = (dVol <= 0.)

                Vol_tol = max(0.5 * GV%Angstrom_H + GV%H_subroundoff, 1e-14 * vol)
                Vol_quit = max(0.9 * GV%Angstrom_H + GV%H_subroundoff, 1e-14 * vol)

                curv_tol = Vol_tol * dV_dL2**2 &
                           * (dV_dL2 * Vol_tol - 2.0 * a * L0 * dVol)
                do_one_L_iter = (a * a * dVol**3) < curv_tol
              endif

              if (use_L0) then
                L(K) = L0
                Vol_err = 0.5*(L(K)*L(K))*(slope + a_3*(3.0-4.0*L(K))) - vol
              elseif (do_one_L_iter) then
                ! One iteration of Newton's method should give an estimate
                ! that is accurate to within Vol_tol.
                L(K) = sqrt(L0*L0 + dVol / dV_dL2)
                Vol_err = 0.5*(L(K)*L(K))*(slope + a_3*(3.0-4.0*L(K))) - vol
              else
                if (dV_dL2*(1.0-L0*L0) < dVol + &
                    dV_dL2 * (Vol_open - Vol)*Ibma_2) then
                  L_max = sqrt(1.0 - (Vol_open - Vol)*Ibma_2)
                else
                  L_max = sqrt(L0*L0 + dVol / dV_dL2)
                endif
                L_min = sqrt(L0*L0 + dVol / (0.5*(slope+a) - a*L_max))

                Vol_err_min = 0.5*(L_min**2)*(slope + a_3*(3.0-4.0*L_min)) - vol
                Vol_err_max = 0.5*(L_max**2)*(slope + a_3*(3.0-4.0*L_max)) - vol
           !    if ((abs(Vol_err_min) <= Vol_quit) .or. (Vol_err_min >= Vol_err_max)) then
                if (abs(Vol_err_min) <= Vol_quit) then
                  L(K) = L_min ; Vol_err = Vol_err_min
                else
                  L(K) = sqrt((L_min**2*Vol_err_max - L_max**2*Vol_err_min) / &
                              (Vol_err_max - Vol_err_min))
                  do itt=1,maxitt
                    Vol_err = 0.5*(L(K)*L(K))*(slope + a_3*(3.0-4.0*L(K))) - vol
                    if (abs(Vol_err) <= Vol_quit) exit
                    ! Take a Newton's method iteration. This equation has proven
                    ! robust enough not to need bracketing.
                    L(K) = L(K) - Vol_err / (L(K)* (slope + a - 2.0*a*L(K)))
                    ! This would be a Newton's method iteration for L^2:
                    !   L(K) = sqrt(L(K)*L(K) - Vol_err / (0.5*(slope+a) - a*L(K)))
                  enddo
                endif ! end of iterative solver
              endif ! end of 1-boundary alternatives.
            endif ! end of a<0 cases.
          endif

          ! Determine the drag contributing to the bottom boundary layer
          ! and the Raleigh drag that acts on each layer.
          if (L(K) > L(K+1)) then
            if (vol_below < bbl_thick) then
              BBL_frac = (1.0-vol_below/bbl_thick)**2
              BBL_visc_frac = BBL_visc_frac + BBL_frac*(L(K) - L(K+1))
            else
              BBL_frac = 0.0
            endif

            if (m==1) then ; Cell_width = G%dy_Cu(I,j)
            else ; Cell_width = G%dx_Cv(i,J) ; endif
            gam = 1.0 - L(K+1)/L(K)
            Rayleigh = US%L_to_Z * CS%cdrag * (L(K)-L(K+1)) * (1.0-BBL_frac) * &
                (12.0*CS%c_Smag*h_vel_pos) /  (12.0*CS%c_Smag*h_vel_pos + &
                 US%L_to_Z*GV%Z_to_H * CS%cdrag * gam*(1.0-gam)*(1.0-1.5*gam) * L(K)**2 * Cell_width)
          else ! This layer feels no drag.
            Rayleigh = 0.0
          endif

          if (m==1) then
            if (Rayleigh > 0.0) then
              v_at_u = set_v_at_u(v, h, G, GV, i, j, k, mask_v, OBC)
              visc%Ray_u(I,j,k) = Rayleigh*sqrt(u(I,j,k)*u(I,j,k) + &
                                                v_at_u*v_at_u + U_bg_sq)
            else ; visc%Ray_u(I,j,k) = 0.0 ; endif
          else
            if (Rayleigh > 0.0) then
              u_at_v = set_u_at_v(u, h, G, GV, i, j, k, mask_u, OBC)
              visc%Ray_v(i,J,k) = Rayleigh*sqrt(v(i,J,k)*v(i,J,k) + &
                                                u_at_v*u_at_v + U_bg_sq)
            else ; visc%Ray_v(i,J,k) = 0.0 ; endif
          endif

        enddo ! k loop to determine L(K).

        ! Set the near-bottom viscosity to a value which will give
        ! the correct stress when the shear occurs over bbl_thick.
        ! See next block for explanation.
        bbl_thick_Z = bbl_thick * GV%H_to_Z
        if (CS%correct_BBL_bounds .and. &
            cdrag_sqrt*ustar(i)*bbl_thick_Z*BBL_visc_frac <= CS%Kv_BBL_min) then
          ! If the bottom stress implies less viscosity than Kv_BBL_min then
          ! set kv_bbl to the bound and recompute bbl_thick to be consistent
          ! but with a ridiculously large upper bound on thickness (for Cd u*=0)
          kv_bbl = CS%Kv_BBL_min
          if (cdrag_sqrt*ustar(i)*BBL_visc_frac*G%Rad_Earth*US%m_to_Z > kv_bbl) then
            bbl_thick_Z = kv_bbl / ( cdrag_sqrt*ustar(i)*BBL_visc_frac )
          else
            bbl_thick_Z = G%Rad_Earth * US%m_to_Z
          endif
        else
          kv_bbl = cdrag_sqrt*ustar(i)*bbl_thick_Z*BBL_visc_frac
        endif

      else ! Not Channel_drag.
        ! Set the near-bottom viscosity to a value which will give
        ! the correct stress when the shear occurs over bbl_thick.
        ! - The bottom stress is tau_b = Cdrag * u_bbl^2
        ! - u_bbl was calculated by averaging flow over CS%Hbbl
        !   (and includes unresolved tidal components)
        ! - u_bbl is embedded in u* since u*^2 = Cdrag u_bbl^2
        ! - The average shear in the BBL is du/dz = 2 * u_bbl / h_bbl
        !   (which assumes a linear profile, hence the "2")
        ! - bbl_thick was bounded to <= 0.5 * CS%Hbbl
        ! - The viscous stress kv_bbl du/dz should balance tau_b
        !      Cdrag u_bbl^2 = kv_bbl du/dz
        !                    = 2 kv_bbl u_bbl
        ! so
        !      kv_bbl = 0.5 h_bbl Cdrag u_bbl
        !             = 0.5 h_bbl sqrt(Cdrag) u*
        bbl_thick_Z = bbl_thick * GV%H_to_Z
        if (CS%correct_BBL_bounds .and. &
            cdrag_sqrt*ustar(i)*bbl_thick_Z <= CS%Kv_BBL_min) then
          ! If the bottom stress implies less viscosity than Kv_BBL_min then
          ! set kv_bbl to the bound and recompute bbl_thick to be consistent
          ! but with a ridiculously large upper bound on thickness (for Cd u*=0)
          kv_bbl = CS%Kv_BBL_min
          if (cdrag_sqrt*ustar(i)*G%Rad_Earth*US%m_to_Z > kv_bbl) then
            bbl_thick_Z = kv_bbl / ( cdrag_sqrt*ustar(i) )
          else
            bbl_thick_Z = G%Rad_Earth * US%m_to_Z
          endif
        else
          kv_bbl = cdrag_sqrt*ustar(i)*bbl_thick_Z
        endif
      endif
      kv_bbl = max(CS%Kv_BBL_min, kv_bbl)
      if (m==1) then
        visc%Kv_bbl_u(I,j) = kv_bbl
        visc%bbl_thick_u(I,j) = bbl_thick_Z
      else
        visc%Kv_bbl_v(i,J) = kv_bbl
        visc%bbl_thick_v(i,J) = bbl_thick_Z
      endif
    endif ; enddo ! end of i loop
  enddo ; enddo ! end of m & j loops

! Offer diagnostics for averaging
  if (CS%id_bbl_thick_u > 0) &
    call post_data(CS%id_bbl_thick_u, visc%bbl_thick_u, CS%diag)
  if (CS%id_kv_bbl_u > 0) &
    call post_data(CS%id_kv_bbl_u, visc%kv_bbl_u, CS%diag)
  if (CS%id_bbl_u > 0) &
    call post_data(CS%id_bbl_u, CS%bbl_u, CS%diag)
  if (CS%id_bbl_thick_v > 0) &
    call post_data(CS%id_bbl_thick_v, visc%bbl_thick_v, CS%diag)
  if (CS%id_kv_bbl_v > 0) &
    call post_data(CS%id_kv_bbl_v, visc%kv_bbl_v, CS%diag)
  if (CS%id_bbl_v > 0) &
    call post_data(CS%id_bbl_v, CS%bbl_v, CS%diag)
  if (CS%id_Ray_u > 0) &
    call post_data(CS%id_Ray_u, visc%Ray_u, CS%diag)
  if (CS%id_Ray_v > 0) &
    call post_data(CS%id_Ray_v, visc%Ray_v, CS%diag)

  if (CS%debug) then
    if (associated(visc%Ray_u) .and. associated(visc%Ray_v)) &
        call uvchksum("Ray [uv]", visc%Ray_u, visc%Ray_v, G%HI, haloshift=0, scale=US%Z_to_m*US%s_to_T)
    if (associated(visc%kv_bbl_u) .and. associated(visc%kv_bbl_v)) &
        call uvchksum("kv_bbl_[uv]", visc%kv_bbl_u, visc%kv_bbl_v, G%HI, &
                      haloshift=0, scale=US%Z2_T_to_m2_s, scalar_pair=.true.)
    if (associated(visc%bbl_thick_u) .and. associated(visc%bbl_thick_v)) &
        call uvchksum("bbl_thick_[uv]", visc%bbl_thick_u, visc%bbl_thick_v, &
                      G%HI, haloshift=0, scale=US%Z_to_m, scalar_pair=.true.)
  endif

end subroutine set_viscous_BBL

!> This subroutine finds a thickness-weighted value of v at the u-points.
function set_v_at_u(v, h, G, GV, i, j, k, mask2dCv, OBC)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in) :: v    !< The meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  integer,                 intent(in) :: i    !< The i-index of the u-location to work on.
  integer,                 intent(in) :: j    !< The j-index of the u-location to work on.
  integer,                 intent(in) :: k    !< The k-index of the u-location to work on.
  real, dimension(SZI_(G),SZJB_(G)),&
                           intent(in) :: mask2dCv !< A multiplicative mask of the v-points
  type(ocean_OBC_type),    pointer    :: OBC  !< A pointer to an open boundary condition structure
  real                                :: set_v_at_u !< The return value of v at u points points in the
                                              !! same units as u, i.e. [L T-1 ~> m s-1] or other units.

  ! This subroutine finds a thickness-weighted value of v at the u-points.
  real :: hwt(0:1,-1:0)    ! Masked weights used to average u onto v [H ~> m or kg m-2].
  real :: hwt_tot          ! The sum of the masked thicknesses [H ~> m or kg m-2].
  integer :: i0, j0, i1, j1

  do j0 = -1,0 ; do i0 = 0,1 ; i1 = i+i0 ; J1 = J+j0
    hwt(i0,j0) = (h(i1,j1,k) + h(i1,j1+1,k)) * mask2dCv(i1,J1)
  enddo ; enddo

  if (associated(OBC)) then ; if (OBC%number_of_segments > 0) then
    do j0 = -1,0 ; do i0 = 0,1 ; if ((OBC%segnum_v(i+i0,J+j0) /= OBC_NONE)) then
      i1 = i+i0 ; J1 = J+j0
      if (OBC%segment(OBC%segnum_v(i1,j1))%direction == OBC_DIRECTION_N) then
        hwt(i0,j0) = 2.0 * h(i1,j1,k) * mask2dCv(i1,J1)
      elseif (OBC%segment(OBC%segnum_v(i1,J1))%direction == OBC_DIRECTION_S) then
        hwt(i0,j0) = 2.0 * h(i1,J1+1,k) * mask2dCv(i1,J1)
      endif
    endif ; enddo ; enddo
  endif ; endif

  hwt_tot = (hwt(0,-1) + hwt(1,0)) + (hwt(1,-1) + hwt(0,0))
  set_v_at_u = 0.0
  if (hwt_tot > 0.0) set_v_at_u = &
          ((hwt(0,0) * v(i,J,k) + hwt(1,-1) * v(i+1,J-1,k)) + &
           (hwt(1,0) * v(i+1,J,k) + hwt(0,-1) * v(i,J-1,k))) / hwt_tot

end function set_v_at_u

!> This subroutine finds a thickness-weighted value of u at the v-points.
function set_u_at_v(u, h, G, GV, i, j, k, mask2dCu, OBC)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: u    !< The zonal velocity [L T-1 ~> m s-1] or other units.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  integer,                 intent(in) :: i    !< The i-index of the u-location to work on.
  integer,                 intent(in) :: j    !< The j-index of the u-location to work on.
  integer,                 intent(in) :: k    !< The k-index of the u-location to work on.
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in) :: mask2dCu !< A multiplicative mask of the u-points
  type(ocean_OBC_type),    pointer    :: OBC  !< A pointer to an open boundary condition structure
  real                                :: set_u_at_v !< The return value of u at v points in the
                                              !! same units as u, i.e. [L T-1 ~> m s-1] or other units.

  ! This subroutine finds a thickness-weighted value of u at the v-points.
  real :: hwt(-1:0,0:1)    ! Masked weights used to average u onto v [H ~> m or kg m-2].
  real :: hwt_tot          ! The sum of the masked thicknesses [H ~> m or kg m-2].
  integer :: i0, j0, i1, j1

  do j0 = 0,1 ; do i0 = -1,0 ; I1 = I+i0 ; j1 = j+j0
    hwt(i0,j0) = (h(i1,j1,k) + h(i1+1,j1,k)) * mask2dCu(I1,j1)
  enddo ; enddo

  if (associated(OBC)) then ; if (OBC%number_of_segments > 0) then
    do j0 = 0,1 ; do i0 = -1,0 ; if ((OBC%segnum_u(I+i0,j+j0) /= OBC_NONE)) then
      I1 = I+i0 ; j1 = j+j0
      if (OBC%segment(OBC%segnum_u(I1,j1))%direction == OBC_DIRECTION_E) then
        hwt(i0,j0) = 2.0 * h(I1,j1,k) * mask2dCu(I1,j1)
      elseif (OBC%segment(OBC%segnum_u(I1,j1))%direction == OBC_DIRECTION_W) then
        hwt(i0,j0) = 2.0 * h(I1+1,j1,k) * mask2dCu(I1,j1)
      endif
    endif ; enddo ; enddo
  endif ; endif

  hwt_tot = (hwt(-1,0) + hwt(0,1)) + (hwt(0,0) + hwt(-1,1))
  set_u_at_v = 0.0
  if (hwt_tot > 0.0) set_u_at_v = &
          ((hwt(0,0) * u(I,j,k) + hwt(-1,1) * u(I-1,j+1,k)) + &
           (hwt(-1,0) * u(I-1,j,k) + hwt(0,1) * u(I,j+1,k))) / hwt_tot

end function set_u_at_v

!> Calculates the thickness of the surface boundary layer for applying an elevated viscosity.
!!
!! A bulk Richardson criterion or the thickness of the topmost NKML layers (with a bulk mixed layer)
!! are currently used.  The thicknesses are given in terms of fractional layers, so that this
!! thickness will move as the thickness of the topmost layers change.
subroutine set_viscous_ML(u, v, h, tv, forces, visc, dt, G, GV, US, CS, symmetrize)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure containing pointers to any available
                                                 !! thermodynamic fields. Absent fields have
                                                 !! NULL ptrs.
  type(mech_forcing),      intent(in)    :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),     intent(inout) :: visc !< A structure containing vertical viscosities and
                                                 !! related fields.
  real,                    intent(in)    :: dt   !< Time increment [T ~> s].
  type(set_visc_CS),       pointer       :: CS   !< The control structure returned by a previous
                                                 !! call to set_visc_init.
  logical,        optional, intent(in)    :: symmetrize !< If present and true, do extra calculations
                                                 !! of those values in visc that would be
                                                 !! calculated with symmetric memory.

  ! Local variables
  real, dimension(SZIB_(G)) :: &
    htot, &     !   The total depth of the layers being that are within the
                ! surface mixed layer [H ~> m or kg m-2].
    Thtot, &    !   The integrated temperature of layers that are within the
                ! surface mixed layer [H degC ~> m degC or kg degC m-2].
    Shtot, &    !   The integrated salt of layers that are within the
                ! surface mixed layer [H ppt ~> m ppt or kg ppt m-2].
    Rhtot, &    !   The integrated density of layers that are within the surface mixed layer
                ! [H R ~> kg m-2 or kg2 m-5].  Rhtot is only used if no
                ! equation of state is used.
    uhtot, &    !   The depth integrated zonal and meridional velocities within
    vhtot, &    ! the surface mixed layer [H L T-1 ~> m2 s-1 or kg m-1 s-1].
    Idecay_len_TKE, & ! The inverse of a turbulence decay length scale [H-1 ~> m-1 or m2 kg-1].
    dR_dT, &    !   Partial derivative of the density at the base of layer nkml
                ! (roughly the base of the mixed layer) with temperature [R degC-1 ~> kg m-3 degC-1].
    dR_dS, &    !   Partial derivative of the density at the base of layer nkml
                ! (roughly the base of the mixed layer) with salinity [R ppt-1 ~> kg m-3 ppt-1].
    ustar, &    !   The surface friction velocity under ice shelves [Z T-1 ~> m s-1].
    press, &    ! The pressure at which dR_dT and dR_dS are evaluated [R L2 T-2 ~> Pa].
    T_EOS, &    ! The potential temperature at which dR_dT and dR_dS are evaluated [degC]
    S_EOS       ! The salinity at which dR_dT and dR_dS are evaluated [ppt].
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    mask_u      ! A mask that disables any contributions from u points that
                ! are land or past open boundary conditions [nondim], 0 or 1.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    mask_v      ! A mask that disables any contributions from v points that
                ! are land or past open boundary conditions [nondim], 0 or 1.
  real :: h_at_vel(SZIB_(G),SZK_(GV))! Layer thickness at velocity points,
                ! using an upwind-biased second order accurate estimate based
                ! on the previous velocity direction [H ~> m or kg m-2].
  integer :: k_massive(SZIB_(G)) ! The k-index of the deepest layer yet found
                ! that has more than h_tiny thickness and will be in the
                ! viscous mixed layer.
  real :: Uh2   ! The squared magnitude of the difference between the velocity
                ! integrated through the mixed layer and the velocity of the
                ! interior layer layer times the depth of the the mixed layer
                ! [H2 Z2 T-2 ~> m4 s-2 or kg2 m-2 s-2].
  real :: htot_vel  ! Sum of the layer thicknesses up to some point [H ~> m or kg m-2].
  real :: hwtot     ! Sum of the thicknesses used to calculate
                    ! the near-bottom velocity magnitude [H ~> m or kg m-2].
  real :: hutot     ! Running sum of thicknesses times the
                    ! velocity magnitudes [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real :: hweight   ! The thickness of a layer that is within Hbbl
                    ! of the bottom [H ~> m or kg m-2].
  real :: tbl_thick_Z  ! The thickness of the top boundary layer [Z ~> m].

  real :: hlay      ! The layer thickness at velocity points [H ~> m or kg m-2].
  real :: I_2hlay   ! 1 / 2*hlay [H-1 ~> m-1 or m2 kg-1].
  real :: T_lay     ! The layer temperature at velocity points [degC].
  real :: S_lay     ! The layer salinity at velocity points [ppt].
  real :: Rlay      ! The layer potential density at velocity points [R ~> kg m-3].
  real :: Rlb       ! The potential density of the layer below [R ~> kg m-3].
  real :: v_at_u    ! The meridonal velocity at a zonal velocity point [L T-1 ~> m s-1].
  real :: u_at_v    ! The zonal velocity at a meridonal velocity point [L T-1 ~> m s-1].
  real :: gHprime   ! The mixed-layer internal gravity wave speed squared, based
                    ! on the mixed layer thickness and density difference across
                    ! the base of the mixed layer [L2 T-2 ~> m2 s-2].
  real :: RiBulk    ! The bulk Richardson number below which water is in the
                    ! viscous mixed layer, including reduction for turbulent
                    ! decay. Nondimensional.
  real :: dt_Rho0   ! The time step divided by the conversion from the layer
                    ! thickness to layer mass [T H Z-1 R-1 ~> s m3 kg-1 or s].
  real :: g_H_Rho0  !   The gravitational acceleration times the conversion from H to m divided
                    ! by the mean density [L2 T-2 H-1 R-1 ~> m4 s-2 kg-1 or m7 s-2 kg-2].
  real :: ustarsq     ! 400 times the square of ustar, times
                      ! Rho0 divided by G_Earth and the conversion
                      ! from m to thickness units [H R ~> kg m-2 or kg2 m-5].
  real :: cdrag_sqrt_Z  ! Square root of the drag coefficient, times a unit conversion
                      ! factor from lateral lengths to vertical depths [Z L-1 ~> 1].
  real :: cdrag_sqrt  ! Square root of the drag coefficient [nondim].
  real :: oldfn       ! The integrated energy required to
                      ! entrain up to the bottom of the layer,
                      ! divided by G_Earth [H R ~> kg m-2 or kg2 m-5].
  real :: Dfn         ! The increment in oldfn for entraining
                      ! the layer [H R ~> kg m-2 or kg2 m-5].
  real :: Dh          ! The increment in layer thickness from
                      ! the present layer [H ~> m or kg m-2].
  real :: U_bg_sq   ! The square of an assumed background velocity, for
                    ! calculating the mean magnitude near the top for use in
                    ! the quadratic surface drag [L2 T-2 ~> m2 s-2].
  real :: h_tiny    ! A very small thickness [H ~> m or kg m-2]. Layers that are less than
                    ! h_tiny can not be the deepest in the viscous mixed layer.
  real :: absf      ! The absolute value of f averaged to velocity points [T-1 ~> s-1].
  real :: U_star    ! The friction velocity at velocity points [Z T-1 ~> m s-1].
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: Rho0x400_G ! 400*Rho0/G_Earth, times unit conversion factors
                     ! [R T2 H Z-2 ~> kg s2 m-4 or kg2 s2 m-7].
                     ! The 400 is a constant proposed by Killworth and Edwards, 1999.
  real :: ustar1    ! ustar [H T-1 ~> m s-1 or kg m-2 s-1]
  real :: h2f2      ! (h*2*f)^2 [H2 T-2 ~> m2 s-2 or kg2 m-4 s-2]
  logical :: use_EOS, do_any, do_any_shelf, do_i(SZIB_(G))
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, K2, nkmb, nkml, n
  type(ocean_OBC_type), pointer :: OBC => NULL()

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nkmb = GV%nk_rho_varies ; nkml = GV%nkml

  if (.not.associated(CS)) call MOM_error(FATAL,"MOM_set_viscosity(visc_ML): "//&
         "Module must be initialized before it is used.")
  if (.not.(CS%dynamic_viscous_ML .or. associated(forces%frac_shelf_u) .or. &
            associated(forces%frac_shelf_v)) ) return

  if (present(symmetrize)) then ; if (symmetrize) then
    Jsq = js-1 ; Isq = is-1
  endif ; endif

  Rho0x400_G = 400.0*(GV%Rho0/(US%L_to_Z**2 * GV%g_Earth)) * GV%Z_to_H
  U_bg_sq = CS%drag_bg_vel * CS%drag_bg_vel
  cdrag_sqrt = sqrt(CS%cdrag)
  cdrag_sqrt_Z = US%L_to_Z * sqrt(CS%cdrag)

  OBC => CS%OBC
  use_EOS = associated(tv%eqn_of_state)
  dt_Rho0 = dt / GV%H_to_RZ
  h_neglect = GV%H_subroundoff
  h_tiny = 2.0*GV%Angstrom_H + h_neglect
  g_H_Rho0 = (GV%g_Earth*GV%H_to_Z) / (GV%Rho0)

  if (associated(forces%frac_shelf_u) .neqv. associated(forces%frac_shelf_v)) &
    call MOM_error(FATAL, "set_viscous_ML: one of forces%frac_shelf_u and "//&
                   "forces%frac_shelf_v is associated, but the other is not.")

  if (associated(forces%frac_shelf_u)) then
    ! This configuration has ice shelves, and the appropriate variables need to be
    ! allocated.  If the arrays have already been allocated, these calls do nothing.
    call safe_alloc_ptr(visc%tauy_shelf, G%isd, G%ied, G%JsdB, G%JedB)
    call safe_alloc_ptr(visc%tbl_thick_shelf_u, G%IsdB, G%IedB, G%jsd, G%jed)
    call safe_alloc_ptr(visc%tbl_thick_shelf_v, G%isd, G%ied, G%JsdB, G%JedB)
    call safe_alloc_ptr(visc%kv_tbl_shelf_u, G%IsdB, G%IedB, G%jsd, G%jed)
    call safe_alloc_ptr(visc%kv_tbl_shelf_v, G%isd, G%ied, G%JsdB, G%JedB)
    call safe_alloc_ptr(visc%taux_shelf, G%IsdB, G%IedB, G%jsd, G%jed)
    call safe_alloc_ptr(visc%tauy_shelf, G%isd, G%ied, G%JsdB, G%JedB)

    !  With a linear drag law under shelves, the friction velocity is already known.
!    if (CS%linear_drag) ustar(:) = cdrag_sqrt_Z*CS%drag_bg_vel
  endif

  !$OMP parallel do default(shared)
  do J=js-1,je ; do i=is-1,ie+1
    mask_v(i,J) = G%mask2dCv(i,J)
  enddo ; enddo
  !$OMP parallel do default(shared)
  do j=js-1,je+1 ; do I=is-1,ie
    mask_u(I,j) = G%mask2dCu(I,j)
  enddo ; enddo

  if (associated(OBC)) then ; do n=1,OBC%number_of_segments
    ! Now project bottom depths across cell-corner points in the OBCs.  The two
    ! projections have to occur in sequence and can not be combined easily.
    if (.not. OBC%segment(n)%on_pe) cycle
    ! Use a one-sided projection of bottom depths at OBC points.
    I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
    if (OBC%segment(n)%is_N_or_S .and. (J >= js-1) .and. (J <= je)) then
      do I = max(is-1,OBC%segment(n)%HI%IsdB), min(ie,OBC%segment(n)%HI%IedB)
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) mask_u(I,j+1) = 0.0
        if (OBC%segment(n)%direction == OBC_DIRECTION_S) mask_u(I,j) = 0.0
      enddo
    elseif (OBC%segment(n)%is_E_or_W .and. (I >= is-1) .and. (I <= je)) then
      do J = max(js-1,OBC%segment(n)%HI%JsdB), min(je,OBC%segment(n)%HI%JedB)
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) mask_v(i+1,J) = 0.0
        if (OBC%segment(n)%direction == OBC_DIRECTION_W) mask_v(i,J) = 0.0
      enddo
    endif
  enddo ; endif

  !$OMP parallel do default(private) shared(u,v,h,tv,forces,visc,dt,G,GV,US,CS,use_EOS,dt_Rho0, &
  !$OMP                                     h_neglect,h_tiny,g_H_Rho0,js,je,OBC,Isq,Ieq,nz,  &
  !$OMP                                     U_bg_sq,mask_v,cdrag_sqrt,cdrag_sqrt_Z,Rho0x400_G,nkml)
  do j=js,je  ! u-point loop
    if (CS%dynamic_viscous_ML) then
      do_any = .false.
      do I=Isq,Ieq
        htot(I) = 0.0
        if (G%mask2dCu(I,j) < 0.5) then
          do_i(I) = .false. ; visc%nkml_visc_u(I,j) = nkml
        else
          do_i(I) = .true. ; do_any = .true.
          k_massive(I) = nkml
          Thtot(I) = 0.0 ; Shtot(I) = 0.0 ; Rhtot(i) = 0.0
          uhtot(I) = dt_Rho0 * forces%taux(I,j)
          vhtot(I) = 0.25 * dt_Rho0 * ((forces%tauy(i,J) + forces%tauy(i+1,J-1)) + &
                                       (forces%tauy(i,J-1) + forces%tauy(i+1,J)))

          if (CS%omega_frac >= 1.0) then ; absf = 2.0*CS%omega ; else
            absf = 0.5*(abs(G%CoriolisBu(I,J)) + abs(G%CoriolisBu(I,J-1)))
            if (CS%omega_frac > 0.0) &
              absf = sqrt(CS%omega_frac*4.0*CS%omega**2 + (1.0-CS%omega_frac)*absf**2)
          endif
          U_star = max(CS%ustar_min, 0.5 * (forces%ustar(i,j) + forces%ustar(i+1,j)))
          Idecay_len_TKE(I) = ((absf / U_star) * CS%TKE_decay) * GV%H_to_Z
        endif
      enddo

      if (do_any) then ; do k=1,nz
        if (k > nkml) then
          do_any = .false.
          if (use_EOS .and. (k==nkml+1)) then
            ! Find dRho/dT and dRho_dS.
            do I=Isq,Ieq
              press(I) = (GV%H_to_RZ*GV%g_Earth) * htot(I)
              if (associated(tv%p_surf)) press(I) = press(I) + 0.5*(tv%p_surf(i,j)+tv%p_surf(i+1,j))
              k2 = max(1,nkml)
              I_2hlay = 1.0 / (h(i,j,k2) + h(i+1,j,k2) + h_neglect)
              T_EOS(I) = (h(i,j,k2)*tv%T(i,j,k2) + h(i+1,j,k2)*tv%T(i+1,j,k2)) * I_2hlay
              S_EOS(I) = (h(i,j,k2)*tv%S(i,j,k2) + h(i+1,j,k2)*tv%S(i+1,j,k2)) * I_2hlay
            enddo
            call calculate_density_derivs(T_EOS, S_EOS, press, dR_dT, dR_dS, tv%eqn_of_state, &
                                          (/Isq-G%IsdB+1,Ieq-G%IsdB+1/) )
          endif

          do I=Isq,Ieq ; if (do_i(I)) then

            hlay = 0.5*(h(i,j,k) + h(i+1,j,k))
            if (hlay > h_tiny) then ! Only consider non-vanished layers.
              I_2hlay = 1.0 / (h(i,j,k) + h(i+1,j,k))
              v_at_u = 0.5 * (h(i,j,k)   * (v(i,J,k) + v(i,J-1,k)) + &
                              h(i+1,j,k) * (v(i+1,J,k) + v(i+1,J-1,k))) * I_2hlay
              Uh2 = ((uhtot(I) - htot(I)*u(I,j,k))**2 + (vhtot(I) - htot(I)*v_at_u)**2)

              if (use_EOS) then
                T_lay = (h(i,j,k)*tv%T(i,j,k) + h(i+1,j,k)*tv%T(i+1,j,k)) * I_2hlay
                S_lay = (h(i,j,k)*tv%S(i,j,k) + h(i+1,j,k)*tv%S(i+1,j,k)) * I_2hlay
                gHprime = g_H_Rho0 * (dR_dT(I) * (T_lay*htot(I) - Thtot(I)) + &
                                      dR_dS(I) * (S_lay*htot(I) - Shtot(I)))
              else
                gHprime = g_H_Rho0 * (GV%Rlay(k)*htot(I) - Rhtot(I))
              endif

              if (gHprime > 0.0) then
                RiBulk = CS%bulk_Ri_ML * exp(-htot(I) * Idecay_len_TKE(I))
                if (RiBulk * Uh2 <= (htot(I)**2) * gHprime) then
                  visc%nkml_visc_u(I,j) = real(k_massive(I))
                  do_i(I) = .false.
                elseif (RiBulk * Uh2 <= (htot(I) + hlay)**2 * gHprime) then
                  visc%nkml_visc_u(I,j) = real(k-1) + &
                    ( sqrt(RiBulk * Uh2 / gHprime) - htot(I) ) / hlay
                  do_i(I) = .false.
                endif
              endif
              k_massive(I) = k
            endif ! hlay > h_tiny

            if (do_i(I)) do_any = .true.
          endif ; enddo

          if (.not.do_any) exit ! All columns are done.
        endif

        do I=Isq,Ieq ; if (do_i(I)) then
          htot(I) = htot(I) + 0.5 * (h(i,j,k) + h(i+1,j,k))
          uhtot(I) = uhtot(I) + 0.5 * (h(i,j,k) + h(i+1,j,k)) * u(I,j,k)
          vhtot(I) = vhtot(I) + 0.25 * (h(i,j,k) * (v(i,J,k) + v(i,J-1,k)) + &
                                        h(i+1,j,k) * (v(i+1,J,k) + v(i+1,J-1,k)))
          if (use_EOS) then
            Thtot(I) = Thtot(I) + 0.5 * (h(i,j,k)*tv%T(i,j,k) + h(i+1,j,k)*tv%T(i+1,j,k))
            Shtot(I) = Shtot(I) + 0.5 * (h(i,j,k)*tv%S(i,j,k) + h(i+1,j,k)*tv%S(i+1,j,k))
          else
            Rhtot(i) = Rhtot(i) + 0.5 * (h(i,j,k) + h(i+1,j,k)) * GV%Rlay(k)
          endif
        endif ; enddo
      enddo ; endif

      if (do_any) then ; do I=Isq,Ieq ; if (do_i(I)) then
        visc%nkml_visc_u(I,j) = k_massive(I)
      endif ; enddo ; endif
    endif ! dynamic_viscous_ML

    do_any_shelf = .false.
    if (associated(forces%frac_shelf_u)) then
      do I=Isq,Ieq
        if (forces%frac_shelf_u(I,j)*G%mask2dCu(I,j) == 0.0) then
          do_i(I) = .false.
          visc%tbl_thick_shelf_u(I,j) = 0.0 ; visc%kv_tbl_shelf_u(I,j) = 0.0
        else
          do_i(I) = .true. ; do_any_shelf = .true.
        endif
      enddo
    endif

    if (do_any_shelf) then
      do k=1,nz ; do I=Isq,Ieq ; if (do_i(I)) then
        if (u(I,j,k) *(h(i+1,j,k) - h(i,j,k)) >= 0) then
          h_at_vel(i,k) = 2.0*h(i,j,k)*h(i+1,j,k) / &
                          (h(i,j,k) + h(i+1,j,k) + h_neglect)
        else
          h_at_vel(i,k) =  0.5 * (h(i,j,k) + h(i+1,j,k))
        endif
      else
        h_at_vel(I,k) = 0.0 ; ustar(I) = 0.0
      endif ; enddo ; enddo

      do I=Isq,Ieq ; if (do_i(I)) then
        htot_vel = 0.0 ; hwtot = 0.0 ; hutot = 0.0
        Thtot(I) = 0.0 ; Shtot(I) = 0.0
        if (use_EOS .or. .not.CS%linear_drag) then ; do k=1,nz
          if (htot_vel>=CS%Htbl_shelf) exit ! terminate the k loop
          hweight = MIN(CS%Htbl_shelf - htot_vel, h_at_vel(i,k))
          if (hweight <= 1.5*GV%Angstrom_H + h_neglect) cycle

          htot_vel  = htot_vel + h_at_vel(i,k)
          hwtot = hwtot + hweight

          if (.not.CS%linear_drag) then
            v_at_u = set_v_at_u(v, h, G, GV, i, j, k, mask_v, OBC)
            hutot = hutot + hweight * sqrt(u(I,j,k)**2 + &
                                           v_at_u**2 + U_bg_sq)
          endif
          if (use_EOS) then
            Thtot(I) = Thtot(I) + hweight * 0.5 * (tv%T(i,j,k) + tv%T(i+1,j,k))
            Shtot(I) = Shtot(I) + hweight * 0.5 * (tv%S(i,j,k) + tv%S(i+1,j,k))
          endif
        enddo ; endif

        if ((.not.CS%linear_drag) .and. (hwtot > 0.0)) then
          ustar(I) = cdrag_sqrt_Z * hutot/hwtot
        else
          ustar(I) = cdrag_sqrt_Z * CS%drag_bg_vel
        endif

        if (use_EOS) then ; if (hwtot > 0.0) then
          T_EOS(I) = Thtot(I)/hwtot ; S_EOS(I) = Shtot(I)/hwtot
        else
          T_EOS(I) = 0.0 ; S_EOS(I) = 0.0
        endif ; endif
      endif ; enddo ! I-loop

      if (use_EOS) then
        call calculate_density_derivs(T_EOS, S_EOS, forces%p_surf(:,j), dR_dT, dR_dS, &
                                      tv%eqn_of_state, (/Isq-G%IsdB+1,Ieq-G%IsdB+1/) )
      endif

      do I=Isq,Ieq ; if (do_i(I)) then
  !  The 400.0 in this expression is the square of a constant proposed
  !  by Killworth and Edwards, 1999, in equation (2.20).
        ustarsq = Rho0x400_G * ustar(i)**2
        htot(i) = 0.0
        if (use_EOS) then
          Thtot(i) = 0.0 ; Shtot(i) = 0.0
          do k=1,nz-1
            if (h_at_vel(i,k) <= 0.0) cycle
            T_Lay = 0.5 * (tv%T(i,j,k) + tv%T(i+1,j,k))
            S_Lay = 0.5 * (tv%S(i,j,k) + tv%S(i+1,j,k))
            oldfn = dR_dT(i)*(T_Lay*htot(i) - Thtot(i)) + dR_dS(i)*(S_Lay*htot(i) - Shtot(i))
            if (oldfn >= ustarsq) exit

            Dfn = (dR_dT(i)*(0.5*(tv%T(i,j,k+1)+tv%T(i+1,j,k+1)) - T_Lay) + &
                   dR_dS(i)*(0.5*(tv%S(i,j,k+1)+tv%S(i+1,j,k+1)) - S_Lay)) * &
                  (h_at_vel(i,k)+htot(i))
            if ((oldfn + Dfn) <= ustarsq) then
              Dh = h_at_vel(i,k)
            else
              Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn) / (Dfn))
            endif

            htot(i) = htot(i) + Dh
            Thtot(i) = Thtot(i) + T_Lay*Dh ; Shtot(i) = Shtot(i) + S_Lay*Dh
          enddo
          if ((oldfn < ustarsq) .and. (h_at_vel(i,nz) > 0.0)) then
            T_Lay = 0.5*(tv%T(i,j,nz) + tv%T(i+1,j,nz))
            S_Lay = 0.5*(tv%S(i,j,nz) + tv%S(i+1,j,nz))
            if (dR_dT(i)*(T_Lay*htot(i) - Thtot(i)) + &
                dR_dS(i)*(S_Lay*htot(i) - Shtot(i)) < ustarsq) &
              htot(i) = htot(i) + h_at_vel(i,nz)
          endif ! Examination of layer nz.
        else  ! Use Rlay as the density variable.
          Rhtot = 0.0
          do k=1,nz-1
            Rlay = GV%Rlay(k) ; Rlb = GV%Rlay(k+1)

            oldfn = Rlay*htot(i) - Rhtot(i)
            if (oldfn >= ustarsq) exit

            Dfn = (Rlb - Rlay)*(h_at_vel(i,k)+htot(i))
            if ((oldfn + Dfn) <= ustarsq) then
              Dh = h_at_vel(i,k)
            else
              Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn) / (Dfn))
            endif

            htot(i) = htot(i) + Dh
            Rhtot(i) = Rhtot(i) + Rlay*Dh
          enddo
          if (GV%Rlay(nz)*htot(i) - Rhtot(i) < ustarsq) &
            htot(i) = htot(i) + h_at_vel(i,nz)
        endif ! use_EOS

       !visc%tbl_thick_shelf_u(I,j) = GV%H_to_Z * max(CS%Htbl_shelf_min, &
       !    htot(I) / (0.5 + sqrt(0.25 + &
       !                 (htot(i)*(G%CoriolisBu(I,J-1)+G%CoriolisBu(I,J)))**2 / &
       !                 (ustar(i)*GV%Z_to_H)**2 )) )
        ustar1 = ustar(i)*GV%Z_to_H
        h2f2 = (htot(i)*(G%CoriolisBu(I,J-1)+G%CoriolisBu(I,J)) + h_neglect*CS%omega)**2
        tbl_thick_Z = GV%H_to_Z * max(CS%Htbl_shelf_min, &
            ( htot(I)*ustar1 ) / ( 0.5*ustar1 + sqrt((0.5*ustar1)**2 + h2f2 ) ) )
        visc%tbl_thick_shelf_u(I,j) = tbl_thick_Z
        visc%Kv_tbl_shelf_u(I,j) = max(CS%Kv_TBL_min, cdrag_sqrt*ustar(i)*tbl_thick_Z)
      endif ; enddo ! I-loop
    endif ! do_any_shelf

  enddo ! j-loop at u-points

  !$OMP parallel do default(private) shared(u,v,h,tv,forces,visc,dt,G,GV,US,CS,use_EOS,dt_Rho0, &
  !$OMP                                     h_neglect,h_tiny,g_H_Rho0,is,ie,OBC,Jsq,Jeq,nz, &
  !$OMP                                     U_bg_sq,cdrag_sqrt,cdrag_sqrt_Z,Rho0x400_G,nkml,mask_u)
  do J=Jsq,Jeq  ! v-point loop
    if (CS%dynamic_viscous_ML) then
      do_any = .false.
      do i=is,ie
        htot(i) = 0.0
        if (G%mask2dCv(i,J) < 0.5) then
          do_i(i) = .false. ; visc%nkml_visc_v(i,J) = nkml
        else
          do_i(i) = .true. ; do_any = .true.
          k_massive(i) = nkml
          Thtot(i) = 0.0 ; Shtot(i) = 0.0 ; Rhtot(i) = 0.0
          vhtot(i) = dt_Rho0 * forces%tauy(i,J)
          uhtot(i) = 0.25 * dt_Rho0 * ((forces%taux(I,j) + forces%taux(I-1,j+1)) + &
                                       (forces%taux(I-1,j) + forces%taux(I,j+1)))

         if (CS%omega_frac >= 1.0) then ; absf = 2.0*CS%omega ; else
           absf = 0.5*(abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))
           if (CS%omega_frac > 0.0) &
             absf = sqrt(CS%omega_frac*4.0*CS%omega**2 + (1.0-CS%omega_frac)*absf**2)
         endif

         U_star = max(CS%ustar_min, 0.5 * (forces%ustar(i,j) + forces%ustar(i,j+1)))
         Idecay_len_TKE(i) = ((absf / U_star) * CS%TKE_decay) * GV%H_to_Z

        endif
      enddo

      if (do_any) then ; do k=1,nz
        if (k > nkml) then
          do_any = .false.
          if (use_EOS .and. (k==nkml+1)) then
            ! Find dRho/dT and dRho_dS.
            do i=is,ie
              press(i) = (GV%H_to_RZ * GV%g_Earth) * htot(i)
              if (associated(tv%p_surf)) press(i) = press(i) + 0.5*(tv%p_surf(i,j)+tv%p_surf(i,j+1))
              k2 = max(1,nkml)
              I_2hlay = 1.0 / (h(i,j,k2) + h(i,j+1,k2) + h_neglect)
              T_EOS(i) = (h(i,j,k2)*tv%T(i,j,k2) + h(i,j+1,k2)*tv%T(i,j+1,k2)) * I_2hlay
              S_EOS(i) = (h(i,j,k2)*tv%S(i,j,k2) + h(i,j+1,k2)*tv%S(i,j+1,k2)) * I_2hlay
            enddo
            call calculate_density_derivs(T_EOS, S_EOS, press, dR_dT, dR_dS, &
                                          tv%eqn_of_state, (/is-G%IsdB+1,ie-G%IsdB+1/) )
          endif

          do i=is,ie ; if (do_i(i)) then

            hlay = 0.5*(h(i,j,k) + h(i,j+1,k))
            if (hlay > h_tiny) then ! Only consider non-vanished layers.
              I_2hlay = 1.0 / (h(i,j,k) + h(i,j+1,k))
              u_at_v = 0.5 * (h(i,j,k)   * (u(I-1,j,k)   + u(I,j,k)) + &
                              h(i,j+1,k) * (u(I-1,j+1,k) + u(I,j+1,k))) * I_2hlay
              Uh2 = ((uhtot(I) - htot(I)*u_at_v)**2 + (vhtot(I) - htot(I)*v(i,J,k))**2)

              if (use_EOS) then
                T_lay = (h(i,j,k)*tv%T(i,j,k) + h(i,j+1,k)*tv%T(i,j+1,k)) * I_2hlay
                S_lay = (h(i,j,k)*tv%S(i,j,k) + h(i,j+1,k)*tv%S(i,j+1,k)) * I_2hlay
                gHprime = g_H_Rho0 * (dR_dT(i) * (T_lay*htot(i) - Thtot(i)) + &
                                      dR_dS(i) * (S_lay*htot(i) - Shtot(i)))
              else
                gHprime = g_H_Rho0 * (GV%Rlay(k)*htot(i) - Rhtot(i))
              endif

              if (gHprime > 0.0) then
                RiBulk = CS%bulk_Ri_ML * exp(-htot(i) * Idecay_len_TKE(i))
                if (RiBulk * Uh2 <= htot(i)**2 * gHprime) then
                  visc%nkml_visc_v(i,J) = real(k_massive(i))
                  do_i(i) = .false.
                elseif (RiBulk * Uh2 <= (htot(i) + hlay)**2 * gHprime) then
                  visc%nkml_visc_v(i,J) = real(k-1) + &
                    ( sqrt(RiBulk * Uh2 / gHprime) - htot(i) ) / hlay
                  do_i(i) = .false.
                endif
              endif
              k_massive(i) = k
            endif ! hlay > h_tiny

            if (do_i(i)) do_any = .true.
          endif ; enddo

          if (.not.do_any) exit ! All columns are done.
        endif

        do i=is,ie ; if (do_i(i)) then
          htot(i) = htot(i) + 0.5 * (h(i,J,k) + h(i,j+1,k))
          vhtot(i) = vhtot(i) + 0.5 * (h(i,j,k) + h(i,j+1,k)) * v(i,J,k)
          uhtot(i) = uhtot(i) + 0.25 * (h(i,j,k) * (u(I-1,j,k) + u(I,j,k)) + &
                                        h(i,j+1,k) * (u(I-1,j+1,k) + u(I,j+1,k)))
          if (use_EOS) then
            Thtot(i) = Thtot(i) + 0.5 * (h(i,j,k)*tv%T(i,j,k) + h(i,j+1,k)*tv%T(i,j+1,k))
            Shtot(i) = Shtot(i) + 0.5 * (h(i,j,k)*tv%S(i,j,k) + h(i,j+1,k)*tv%S(i,j+1,k))
          else
            Rhtot(i) = Rhtot(i) + 0.5 * (h(i,j,k) + h(i,j+1,k)) * GV%Rlay(k)
          endif
        endif ; enddo
      enddo ; endif

      if (do_any) then ; do i=is,ie ; if (do_i(i)) then
        visc%nkml_visc_v(i,J) = k_massive(i)
      endif ; enddo ; endif
    endif ! dynamic_viscous_ML

    do_any_shelf = .false.
    if (associated(forces%frac_shelf_v)) then
      do i=is,ie
        if (forces%frac_shelf_v(i,J)*G%mask2dCv(i,J) == 0.0) then
          do_i(i) = .false.
          visc%tbl_thick_shelf_v(i,J) = 0.0 ; visc%kv_tbl_shelf_v(i,J) = 0.0
        else
          do_i(i) = .true. ; do_any_shelf = .true.
        endif
      enddo
    endif

    if (do_any_shelf) then
      do k=1,nz ; do i=is,ie ; if (do_i(i)) then
        if (v(i,J,k) * (h(i,j+1,k) - h(i,j,k)) >= 0) then
          h_at_vel(i,k) = 2.0*h(i,j,k)*h(i,j+1,k) / &
                          (h(i,j,k) + h(i,j+1,k) + h_neglect)
        else
          h_at_vel(i,k) =  0.5 * (h(i,j,k) + h(i,j+1,k))
        endif
      else
        h_at_vel(I,k) = 0.0 ; ustar(i) = 0.0
      endif ; enddo ; enddo

      do i=is,ie ; if (do_i(i)) then
        htot_vel = 0.0 ; hwtot = 0.0 ; hutot = 0.0
        Thtot(i) = 0.0 ; Shtot(i) = 0.0
        if (use_EOS .or. .not.CS%linear_drag) then ; do k=1,nz
          if (htot_vel>=CS%Htbl_shelf) exit ! terminate the k loop
          hweight = MIN(CS%Htbl_shelf - htot_vel, h_at_vel(i,k))
          if (hweight <= 1.5*GV%Angstrom_H + h_neglect) cycle

          htot_vel  = htot_vel + h_at_vel(i,k)
          hwtot = hwtot + hweight

          if (.not.CS%linear_drag) then
            u_at_v = set_u_at_v(u, h, G, GV, i, J, k, mask_u, OBC)
            hutot = hutot + hweight * sqrt(v(i,J,k)**2 + &
                                           u_at_v**2 + U_bg_sq)
          endif
          if (use_EOS) then
            Thtot(i) = Thtot(i) + hweight * 0.5 * (tv%T(i,j,k) + tv%T(i,j+1,k))
            Shtot(i) = Shtot(i) + hweight * 0.5 * (tv%S(i,j,k) + tv%S(i,j+1,k))
          endif
        enddo ; endif

        if (.not.CS%linear_drag) then ; if (hwtot > 0.0) then
          ustar(i) = cdrag_sqrt_Z * hutot/hwtot
        else
          ustar(i) = cdrag_sqrt_Z * CS%drag_bg_vel
        endif ; endif

        if (use_EOS) then ; if (hwtot > 0.0) then
          T_EOS(i) = Thtot(i)/hwtot ; S_EOS(i) = Shtot(i)/hwtot
        else
          T_EOS(i) = 0.0 ; S_EOS(i) = 0.0
        endif ; endif
      endif ; enddo ! I-loop

      if (use_EOS) then
        call calculate_density_derivs(T_EOS, S_EOS, forces%p_surf(:,j), dR_dT, dR_dS, &
                                      tv%eqn_of_state, (/is-G%IsdB+1,ie-G%IsdB+1/) )
      endif

      do i=is,ie ; if (do_i(i)) then
  !  The 400.0 in this expression is the square of a constant proposed
  !  by Killworth and Edwards, 1999, in equation (2.20).
        ustarsq = Rho0x400_G * ustar(i)**2
        htot(i) = 0.0
        if (use_EOS) then
          Thtot(i) = 0.0 ; Shtot(i) = 0.0
          do k=1,nz-1
            if (h_at_vel(i,k) <= 0.0) cycle
            T_Lay = 0.5 * (tv%T(i,j,k) + tv%T(i,j+1,k))
            S_Lay = 0.5 * (tv%S(i,j,k) + tv%S(i,j+1,k))
            oldfn = dR_dT(i)*(T_Lay*htot(i) - Thtot(i)) + dR_dS(i)*(S_Lay*htot(i) - Shtot(i))
            if (oldfn >= ustarsq) exit

            Dfn = (dR_dT(i)*(0.5*(tv%T(i,j,k+1)+tv%T(i,j+1,k+1)) - T_Lay) + &
                   dR_dS(i)*(0.5*(tv%S(i,j,k+1)+tv%S(i,j+1,k+1)) - S_Lay)) * &
                  (h_at_vel(i,k)+htot(i))
            if ((oldfn + Dfn) <= ustarsq) then
              Dh = h_at_vel(i,k)
            else
              Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn) / (Dfn))
            endif

            htot(i) = htot(i) + Dh
            Thtot(i) = Thtot(i) + T_Lay*Dh ; Shtot(i) = Shtot(i) + S_Lay*Dh
          enddo
          if ((oldfn < ustarsq) .and. (h_at_vel(i,nz) > 0.0)) then
            T_Lay = 0.5*(tv%T(i,j,nz) + tv%T(i,j+1,nz))
            S_Lay = 0.5*(tv%S(i,j,nz) + tv%S(i,j+1,nz))
            if (dR_dT(i)*(T_Lay*htot(i) - Thtot(i)) + &
                dR_dS(i)*(S_Lay*htot(i) - Shtot(i)) < ustarsq) &
              htot(i) = htot(i) + h_at_vel(i,nz)
          endif ! Examination of layer nz.
        else  ! Use Rlay as the density variable.
          Rhtot = 0.0
          do k=1,nz-1
            Rlay = GV%Rlay(k) ; Rlb = GV%Rlay(k+1)

            oldfn = Rlay*htot(i) - Rhtot(i)
            if (oldfn >= ustarsq) exit

            Dfn = (Rlb - Rlay)*(h_at_vel(i,k)+htot(i))
            if ((oldfn + Dfn) <= ustarsq) then
              Dh = h_at_vel(i,k)
            else
              Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn) / (Dfn))
            endif

            htot(i) = htot(i) + Dh
            Rhtot = Rhtot + Rlay*Dh
          enddo
          if (GV%Rlay(nz)*htot(i) - Rhtot(i) < ustarsq) &
            htot(i) = htot(i) + h_at_vel(i,nz)
        endif ! use_EOS

       !visc%tbl_thick_shelf_v(i,J) = GV%H_to_Z * max(CS%Htbl_shelf_min, &
       !    htot(i) / (0.5 + sqrt(0.25 + &
       !        (htot(i)*(G%CoriolisBu(I-1,J)+G%CoriolisBu(I,J)))**2 / &
       !        (ustar(i)*GV%Z_to_H)**2 )) )
        ustar1 = ustar(i)*GV%Z_to_H
        h2f2 = (htot(i)*(G%CoriolisBu(I-1,J)+G%CoriolisBu(I,J)) + h_neglect*CS%omega)**2
        tbl_thick_Z = GV%H_to_Z * max(CS%Htbl_shelf_min, &
            ( htot(i)*ustar1 ) / ( 0.5*ustar1 + sqrt((0.5*ustar1)**2 + h2f2 ) ) )
        visc%tbl_thick_shelf_v(i,J) = tbl_thick_Z
        visc%Kv_tbl_shelf_v(i,J) = max(CS%Kv_TBL_min, cdrag_sqrt*ustar(i)*tbl_thick_Z)

      endif ; enddo ! i-loop
    endif ! do_any_shelf

  enddo ! J-loop at v-points

  if (CS%debug) then
    if (associated(visc%nkml_visc_u) .and. associated(visc%nkml_visc_v)) &
      call uvchksum("nkml_visc_[uv]", visc%nkml_visc_u, visc%nkml_visc_v, &
                    G%HI, haloshift=0, scalar_pair=.true.)
  endif
  if (CS%id_nkml_visc_u > 0) &
    call post_data(CS%id_nkml_visc_u, visc%nkml_visc_u, CS%diag)
  if (CS%id_nkml_visc_v > 0) &
    call post_data(CS%id_nkml_visc_v, visc%nkml_visc_v, CS%diag)

end subroutine set_viscous_ML

!> Register any fields associated with the vertvisc_type.
subroutine set_visc_register_restarts(HI, GV, param_file, visc, restart_CS)
  type(hor_index_type),    intent(in)    :: HI         !< A horizontal index type structure.
  type(verticalGrid_type), intent(in)    :: GV         !< The ocean's vertical grid structure.
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                       !! parameters.
  type(vertvisc_type),     intent(inout) :: visc       !< A structure containing vertical
                                                       !! viscosities and related fields.
                                                       !! Allocated here.
  type(MOM_restart_CS),    pointer       :: restart_CS !< A pointer to the restart control structure.
  ! Local variables
  logical :: use_kappa_shear, KS_at_vertex
  logical :: adiabatic, useKPP, useEPBL
  logical :: use_CVMix_shear, MLE_use_PBL_MLD, use_CVMix_conv
  integer :: isd, ied, jsd, jed, nz
  real :: hfreeze !< If hfreeze > 0 [m], melt potential will be computed.
  character(len=40)  :: mdl = "MOM_set_visc"  ! This module's name.
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  call get_param(param_file, mdl, "ADIABATIC", adiabatic, default=.false., &
                 do_not_log=.true.)

  use_kappa_shear = .false. ; KS_at_vertex = .false. ; use_CVMix_shear = .false.
  useKPP = .false. ; useEPBL = .false. ; use_CVMix_conv = .false.

  if (.not.adiabatic) then
    use_kappa_shear = kappa_shear_is_used(param_file)
    KS_at_vertex    = kappa_shear_at_vertex(param_file)
    use_CVMix_shear = CVMix_shear_is_used(param_file)
    use_CVMix_conv  = CVMix_conv_is_used(param_file)
    call get_param(param_file, mdl, "USE_KPP", useKPP, &
                 "If true, turns on the [CVMix] KPP scheme of Large et al., 1994, "//&
                 "to calculate diffusivities and non-local transport in the OBL.", &
                 default=.false., do_not_log=.true.)
    call get_param(param_file, mdl, "ENERGETICS_SFC_PBL", useEPBL, &
                 "If true, use an implied energetics planetary boundary "//&
                 "layer scheme to determine the diffusivity and viscosity "//&
                 "in the surface boundary layer.", default=.false., do_not_log=.true.)
  endif

  if (use_kappa_shear .or. useKPP .or. useEPBL .or. use_CVMix_shear .or. use_CVMix_conv) then
    call safe_alloc_ptr(visc%Kd_shear, isd, ied, jsd, jed, nz+1)
    call register_restart_field(visc%Kd_shear, "Kd_shear", .false., restart_CS, &
                  "Shear-driven turbulent diffusivity at interfaces", "m2 s-1", z_grid='i')
  endif
  if (useKPP .or. useEPBL .or. use_CVMix_shear .or. use_CVMix_conv .or. &
      (use_kappa_shear .and. .not.KS_at_vertex )) then
    call safe_alloc_ptr(visc%Kv_shear, isd, ied, jsd, jed, nz+1)
    call register_restart_field(visc%Kv_shear, "Kv_shear", .false., restart_CS, &
                  "Shear-driven turbulent viscosity at interfaces", "m2 s-1", z_grid='i')
  endif
  if (use_kappa_shear .and. KS_at_vertex) then
    call safe_alloc_ptr(visc%TKE_turb, HI%IsdB, HI%IedB, HI%JsdB, HI%JedB, nz+1)
    call safe_alloc_ptr(visc%Kv_shear_Bu, HI%IsdB, HI%IedB, HI%JsdB, HI%JedB, nz+1)
    call register_restart_field(visc%Kv_shear_Bu, "Kv_shear_Bu", .false., restart_CS, &
                  "Shear-driven turbulent viscosity at vertex interfaces", "m2 s-1", &
                  hor_grid="Bu", z_grid='i')
  elseif (use_kappa_shear) then
    call safe_alloc_ptr(visc%TKE_turb, isd, ied, jsd, jed, nz+1)
  endif

  if (useKPP) then
    ! MOM_bkgnd_mixing uses Kv_slow when KPP is defined.
    call safe_alloc_ptr(visc%Kv_slow, isd, ied, jsd, jed, nz+1)
  endif

  ! visc%MLD is used to communicate the state of the (e)PBL or KPP to the rest of the model
  call get_param(param_file, mdl, "MLE_USE_PBL_MLD", MLE_use_PBL_MLD, &
                 default=.false., do_not_log=.true.)
  ! visc%MLD needs to be allocated when melt potential is computed (HFREEZE>0)
  call get_param(param_file, mdl, "HFREEZE", hfreeze, &
                 default=-1.0, do_not_log=.true.)

  if (MLE_use_PBL_MLD) then
    call safe_alloc_ptr(visc%MLD, isd, ied, jsd, jed)
    call register_restart_field(visc%MLD, "MLD", .false., restart_CS, &
                  "Instantaneous active mixing layer depth", "m")
  endif

  if (hfreeze >= 0.0 .and. .not.MLE_use_PBL_MLD) then
    call safe_alloc_ptr(visc%MLD, isd, ied, jsd, jed)
  endif


end subroutine set_visc_register_restarts

!> Initializes the MOM_set_visc control structure
subroutine set_visc_init(Time, G, GV, US, param_file, diag, visc, CS, restart_CS, OBC)
  type(time_type), target, intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(vertvisc_type),     intent(inout) :: visc !< A structure containing vertical viscosities and
                                                 !! related fields.  Allocated here.
  type(set_visc_CS),       pointer       :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module
  type(MOM_restart_CS),    pointer       :: restart_CS !< A pointer to the restart control structure.
  type(ocean_OBC_type),    pointer       :: OBC  !< A pointer to an open boundary condition structure

  ! Local variables
  real    :: Csmag_chan_dflt, smag_const1, TKE_decay_dflt, bulk_Ri_ML_dflt
  real    :: Kv_background
  real    :: omega_frac_dflt
  real    :: Z_rescale     ! A rescaling factor for heights from the representation in
                           ! a restart file to the internal representation in this run.
  real    :: I_T_rescale   ! A rescaling factor for time from the internal representation in this run
                           ! to the representation in a restart file.
  real    :: Z2_T_rescale  ! A rescaling factor for vertical diffusivities and viscosities from the
                           ! representation in a restart file to the internal representation in this run.
  integer :: i, j, k, is, ie, js, je, n
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  logical :: default_2018_answers
  logical :: use_kappa_shear, adiabatic, use_omega, MLE_use_PBL_MLD
  logical :: use_KPP
  logical :: use_regridding
  character(len=200) :: filename, tideamp_file
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to OBC segment type
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_set_visc"  ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "set_visc_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%OBC => OBC

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%diag => diag

  ! Set default, read and log parameters
  call log_version(param_file, mdl, version, "")
  CS%RiNo_mix = .false.
  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, default=".")
  CS%inputdir = slasher(CS%inputdir)
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  call get_param(param_file, mdl, "SET_VISC_2018_ANSWERS", CS%answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
  call get_param(param_file, mdl, "BOTTOMDRAGLAW", CS%bottomdraglaw, &
                 "If true, the bottom stress is calculated with a drag "//&
                 "law of the form c_drag*|u|*u. The velocity magnitude "//&
                 "may be an assumed value or it may be based on the "//&
                 "actual velocity in the bottommost HBBL, depending on "//&
                 "LINEAR_DRAG.", default=.true.)
  call get_param(param_file, mdl, "CHANNEL_DRAG", CS%Channel_drag, &
                 "If true, the bottom drag is exerted directly on each "//&
                 "layer proportional to the fraction of the bottom it "//&
                 "overlies.", default=.false.)
  call get_param(param_file, mdl, "LINEAR_DRAG", CS%linear_drag, &
                 "If LINEAR_DRAG and BOTTOMDRAGLAW are defined the drag "//&
                 "law is cdrag*DRAG_BG_VEL*u.", default=.false.)
  call get_param(param_file, mdl, "ADIABATIC", adiabatic, default=.false., &
                 do_not_log=.true.)
  if (adiabatic) then
    call log_param(param_file, mdl, "ADIABATIC",adiabatic, &
                 "There are no diapycnal mass fluxes if ADIABATIC is "//&
                 "true. This assumes that KD = KDML = 0.0 and that "//&
                 "there is no buoyancy forcing, but makes the model "//&
                 "faster by eliminating subroutine calls.", default=.false.)
  endif

  if (.not.adiabatic) then
    CS%RiNo_mix = kappa_shear_is_used(param_file)
  endif

  call get_param(param_file, mdl, "PRANDTL_TURB", visc%Prandtl_turb, &
                 "The turbulent Prandtl number applied to shear "//&
                 "instability.", units="nondim", default=1.0)
  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false.)

  call get_param(param_file, mdl, "DYNAMIC_VISCOUS_ML", CS%dynamic_viscous_ML, &
                 "If true, use a bulk Richardson number criterion to "//&
                 "determine the mixed layer thickness for viscosity.", &
                 default=.false.)
  if (CS%dynamic_viscous_ML) then
    call get_param(param_file, mdl, "BULK_RI_ML", bulk_Ri_ML_dflt, default=0.0)
    call get_param(param_file, mdl, "BULK_RI_ML_VISC", CS%bulk_Ri_ML, &
                 "The efficiency with which mean kinetic energy released "//&
                 "by mechanically forced entrainment of the mixed layer "//&
                 "is converted to turbulent kinetic energy.  By default, "//&
                 "BULK_RI_ML_VISC = BULK_RI_ML or 0.", units="nondim", &
                 default=bulk_Ri_ML_dflt)
    call get_param(param_file, mdl, "TKE_DECAY", TKE_decay_dflt, default=0.0)
    call get_param(param_file, mdl, "TKE_DECAY_VISC", CS%TKE_decay, &
                 "TKE_DECAY_VISC relates the vertical rate of decay of "//&
                 "the TKE available for mechanical entrainment to the "//&
                 "natural Ekman depth for use in calculating the dynamic "//&
                 "mixed layer viscosity.  By default, "//&
                 "TKE_DECAY_VISC = TKE_DECAY or 0.", units="nondim", &
                 default=TKE_decay_dflt)
    call get_param(param_file, mdl, "ML_USE_OMEGA", use_omega, &
                 "If true, use the absolute rotation rate instead of the "//&
                 "vertical component of rotation when setting the decay "//&
                   "scale for turbulence.", default=.false., do_not_log=.true.)
    omega_frac_dflt = 0.0
    if (use_omega) then
      call MOM_error(WARNING, "ML_USE_OMEGA is depricated; use ML_OMEGA_FRAC=1.0 instead.")
      omega_frac_dflt = 1.0
    endif
    call get_param(param_file, mdl, "ML_OMEGA_FRAC", CS%omega_frac, &
                   "When setting the decay scale for turbulence, use this "//&
                   "fraction of the absolute rotation rate blended with the "//&
                   "local value of f, as sqrt((1-of)*f^2 + of*4*omega^2).", &
                   units="nondim", default=omega_frac_dflt)
    call get_param(param_file, mdl, "OMEGA", CS%omega, &
                 "The rotation rate of the earth.", units="s-1", &
                 default=7.2921e-5, scale=US%T_to_s)
    ! This give a minimum decay scale that is typically much less than Angstrom.
    CS%ustar_min = 2e-4*CS%omega*(GV%Angstrom_Z + GV%H_to_Z*GV%H_subroundoff)
  else
    call get_param(param_file, mdl, "OMEGA", CS%omega, &
                 "The rotation rate of the earth.", units="s-1", &
                 default=7.2921e-5, scale=US%T_to_s)
  endif

  call get_param(param_file, mdl, "HBBL", CS%Hbbl, &
                 "The thickness of a bottom boundary layer with a "//&
                 "viscosity of KVBBL if BOTTOMDRAGLAW is not defined, or "//&
                 "the thickness over which near-bottom velocities are "//&
                 "averaged for the drag law if BOTTOMDRAGLAW is defined "//&
                 "but LINEAR_DRAG is not.", units="m", fail_if_missing=.true.) ! Rescaled later
  if (CS%bottomdraglaw) then
    call get_param(param_file, mdl, "CDRAG", CS%cdrag, &
                 "CDRAG is the drag coefficient relating the magnitude of "//&
                 "the velocity field to the bottom stress. CDRAG is only "//&
                 "used if BOTTOMDRAGLAW is defined.", units="nondim", &
                 default=0.003)
    call get_param(param_file, mdl, "BBL_USE_TIDAL_BG", CS%BBL_use_tidal_bg, &
                 "Flag to use the tidal RMS amplitude in place of constant "//&
                 "background velocity for computing u* in the BBL. "//&
                 "This flag is only used when BOTTOMDRAGLAW is true and "//&
                 "LINEAR_DRAG is false.", default=.false.)
    if (CS%BBL_use_tidal_bg) then
      call get_param(param_file, mdl, "TIDEAMP_FILE", tideamp_file, &
                   "The path to the file containing the spatially varying "//&
                   "tidal amplitudes with INT_TIDE_DISSIPATION.", default="tideamp.nc")
    else
      call get_param(param_file, mdl, "DRAG_BG_VEL", CS%drag_bg_vel, &
                   "DRAG_BG_VEL is either the assumed bottom velocity (with "//&
                   "LINEAR_DRAG) or an unresolved  velocity that is "//&
                   "combined with the resolved velocity to estimate the "//&
                   "velocity magnitude.  DRAG_BG_VEL is only used when "//&
                   "BOTTOMDRAGLAW is defined.", units="m s-1", default=0.0, scale=US%m_s_to_L_T)
    endif
    call get_param(param_file, mdl, "USE_REGRIDDING", use_regridding, &
         do_not_log = .true., default = .false. )
    call get_param(param_file, mdl, "BBL_USE_EOS", CS%BBL_use_EOS, &
                 "If true, use the equation of state in determining the "//&
                 "properties of the bottom boundary layer.  Otherwise use "//&
                 "the layer target potential densities.  The default of "//&
                 "this is determined by USE_REGRIDDING.", default=use_regridding)
    if (use_regridding .and. (.not. CS%BBL_use_EOS)) &
      call MOM_error(FATAL,"When using MOM6 in ALE mode it is required to "//&
           "set BBL_USE_EOS to True")
  endif
  call get_param(param_file, mdl, "BBL_THICK_MIN", CS%BBL_thick_min, &
                 "The minimum bottom boundary layer thickness that can be "//&
                 "used with BOTTOMDRAGLAW. This might be "//&
                 "Kv/(cdrag*drag_bg_vel) to give Kv as the minimum "//&
                 "near-bottom viscosity.", units="m", default=0.0)  ! Rescaled later
  call get_param(param_file, mdl, "HTBL_SHELF_MIN", CS%Htbl_shelf_min, &
                 "The minimum top boundary layer thickness that can be "//&
                 "used with BOTTOMDRAGLAW. This might be "//&
                 "Kv/(cdrag*drag_bg_vel) to give Kv as the minimum "//&
                 "near-top viscosity.", units="m", default=CS%BBL_thick_min, scale=GV%m_to_H)
  call get_param(param_file, mdl, "HTBL_SHELF", CS%Htbl_shelf, &
                 "The thickness over which near-surface velocities are "//&
                 "averaged for the drag law under an ice shelf.  By "//&
                 "default this is the same as HBBL", units="m", default=CS%Hbbl, scale=GV%m_to_H)
  ! These unit conversions are out outside the get_param calls because the are also defaults.
  CS%Hbbl = CS%Hbbl * GV%m_to_H                   ! Rescale
  CS%BBL_thick_min = CS%BBL_thick_min * GV%m_to_H ! Rescale

  call get_param(param_file, mdl, "KV", Kv_background, &
                 "The background kinematic viscosity in the interior. "//&
                 "The molecular value, ~1e-6 m2 s-1, may be used.", &
                 units="m2 s-1", fail_if_missing=.true.)

  call get_param(param_file, mdl, "USE_KPP", use_KPP, &
                 "If true, turns on the [CVMix] KPP scheme of Large et al., 1994, "//&
                 "to calculate diffusivities and non-local transport in the OBL.", &
                 do_not_log=.true., default=.false.)

  call get_param(param_file, mdl, "KV_BBL_MIN", CS%KV_BBL_min, &
                 "The minimum viscosities in the bottom boundary layer.", &
                 units="m2 s-1", default=Kv_background, scale=US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "KV_TBL_MIN", CS%KV_TBL_min, &
                 "The minimum viscosities in the top boundary layer.", &
                 units="m2 s-1", default=Kv_background, scale=US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "CORRECT_BBL_BOUNDS", CS%correct_BBL_bounds, &
                 "If true, uses the correct bounds on the BBL thickness and "//&
                 "viscosity so that the bottom layer feels the intended drag.", &
                 default=.false.)

  if (CS%Channel_drag) then
    call get_param(param_file, mdl, "SMAG_LAP_CONST", smag_const1, default=-1.0)

    cSmag_chan_dflt = 0.15
    if (smag_const1 >= 0.0) cSmag_chan_dflt = smag_const1

    call get_param(param_file, mdl, "SMAG_CONST_CHANNEL", CS%c_Smag, &
                 "The nondimensional Laplacian Smagorinsky constant used "//&
                 "in calculating the channel drag if it is enabled.  The "//&
                 "default is to use the same value as SMAG_LAP_CONST if "//&
                 "it is defined, or 0.15 if it is not. The value used is "//&
                 "also 0.15 if the specified value is negative.", &
                 units="nondim", default=cSmag_chan_dflt)
    if (CS%c_Smag < 0.0) CS%c_Smag = 0.15
  endif

  call get_param(param_file, mdl, "MLE_USE_PBL_MLD", MLE_use_PBL_MLD, &
                 default=.false., do_not_log=.true.)

  if (CS%RiNo_mix .and. kappa_shear_at_vertex(param_file)) then
    ! This is necessary for reproduciblity across restarts in non-symmetric mode.
    call pass_var(visc%Kv_shear_Bu, G%Domain, position=CORNER, complete=.true.)
  endif

  if (CS%bottomdraglaw) then
    allocate(visc%bbl_thick_u(IsdB:IedB,jsd:jed)) ; visc%bbl_thick_u(:,:) = 0.0
    allocate(visc%kv_bbl_u(IsdB:IedB,jsd:jed)) ; visc%kv_bbl_u(:,:) = 0.0
    allocate(visc%bbl_thick_v(isd:ied,JsdB:JedB)) ; visc%bbl_thick_v(:,:) = 0.0
    allocate(visc%kv_bbl_v(isd:ied,JsdB:JedB)) ; visc%kv_bbl_v(:,:) = 0.0
    allocate(visc%ustar_bbl(isd:ied,jsd:jed)) ; visc%ustar_bbl(:,:) = 0.0
    allocate(visc%TKE_bbl(isd:ied,jsd:jed)) ; visc%TKE_bbl(:,:) = 0.0

    CS%id_bbl_thick_u = register_diag_field('ocean_model', 'bbl_thick_u', &
       diag%axesCu1, Time, 'BBL thickness at u points', 'm', conversion=US%Z_to_m)
    CS%id_kv_bbl_u = register_diag_field('ocean_model', 'kv_bbl_u', diag%axesCu1, &
       Time, 'BBL viscosity at u points', 'm2 s-1', conversion=US%Z2_T_to_m2_s)
    CS%id_bbl_u = register_diag_field('ocean_model', 'bbl_u', diag%axesCu1, &
       Time, 'BBL mean u current', 'm s-1', conversion=US%L_T_to_m_s)
    if (CS%id_bbl_u>0) then
      allocate(CS%bbl_u(IsdB:IedB,jsd:jed)) ; CS%bbl_u(:,:) = 0.0
    endif
    CS%id_bbl_thick_v = register_diag_field('ocean_model', 'bbl_thick_v', &
       diag%axesCv1, Time, 'BBL thickness at v points', 'm', conversion=US%Z_to_m)
    CS%id_kv_bbl_v = register_diag_field('ocean_model', 'kv_bbl_v', diag%axesCv1, &
       Time, 'BBL viscosity at v points', 'm2 s-1', conversion=US%Z2_T_to_m2_s)
    CS%id_bbl_v = register_diag_field('ocean_model', 'bbl_v', diag%axesCv1, &
       Time, 'BBL mean v current', 'm s-1', conversion=US%L_T_to_m_s)
    if (CS%id_bbl_v>0) then
      allocate(CS%bbl_v(isd:ied,JsdB:JedB)) ; CS%bbl_v(:,:) = 0.0
    endif
    if (CS%BBL_use_tidal_bg) then
      allocate(CS%tideamp(isd:ied,jsd:jed)) ; CS%tideamp(:,:) = 0.0
      filename = trim(CS%inputdir) // trim(tideamp_file)
      call log_param(param_file, mdl, "INPUTDIR/TIDEAMP_FILE", filename)
      call MOM_read_data(filename, 'tideamp', CS%tideamp, G%domain, scale=US%m_to_Z*US%T_to_s)
      call pass_var(CS%tideamp,G%domain)
    endif
  endif
  if (CS%Channel_drag) then
    allocate(visc%Ray_u(IsdB:IedB,jsd:jed,nz)) ; visc%Ray_u(:,:,:) = 0.0
    allocate(visc%Ray_v(isd:ied,JsdB:JedB,nz)) ; visc%Ray_v(:,:,:) = 0.0
    CS%id_Ray_u = register_diag_field('ocean_model', 'Rayleigh_u', diag%axesCuL, &
       Time, 'Rayleigh drag velocity at u points', 'm s-1', conversion=US%Z_to_m*US%s_to_T)
    CS%id_Ray_v = register_diag_field('ocean_model', 'Rayleigh_v', diag%axesCvL, &
       Time, 'Rayleigh drag velocity at v points', 'm s-1', conversion=US%Z_to_m*US%s_to_T)
  endif


  if (CS%dynamic_viscous_ML) then
    allocate(visc%nkml_visc_u(IsdB:IedB,jsd:jed)) ; visc%nkml_visc_u(:,:) = 0.0
    allocate(visc%nkml_visc_v(isd:ied,JsdB:JedB)) ; visc%nkml_visc_v(:,:) = 0.0
    CS%id_nkml_visc_u = register_diag_field('ocean_model', 'nkml_visc_u', &
       diag%axesCu1, Time, 'Number of layers in viscous mixed layer at u points', 'm')
    CS%id_nkml_visc_v = register_diag_field('ocean_model', 'nkml_visc_v', &
       diag%axesCv1, Time, 'Number of layers in viscous mixed layer at v points', 'm')
  endif

  call register_restart_field_as_obsolete('Kd_turb','Kd_shear', restart_CS)
  call register_restart_field_as_obsolete('Kv_turb','Kv_shear', restart_CS)

  ! Account for possible changes in dimensional scaling for variables that have been
  ! read from a restart file.
  Z_rescale = 1.0
  if ((US%m_to_Z_restart /= 0.0) .and. (US%m_to_Z_restart /= US%m_to_Z)) &
    Z_rescale = US%m_to_Z / US%m_to_Z_restart
  I_T_rescale = 1.0
  if ((US%s_to_T_restart /= 0.0) .and. (US%s_to_T_restart /= US%s_to_T)) &
    I_T_rescale = US%s_to_T_restart / US%s_to_T
  Z2_T_rescale = Z_rescale**2*I_T_rescale

  if (Z2_T_rescale /= 1.0) then
    if (associated(visc%Kd_shear)) then ; if (query_initialized(visc%Kd_shear, "Kd_shear", restart_CS)) then
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        visc%Kd_shear(i,j,k) = Z2_T_rescale * visc%Kd_shear(i,j,k)
      enddo ; enddo ; enddo
    endif ; endif

    if (associated(visc%Kv_shear)) then ; if (query_initialized(visc%Kv_shear, "Kv_shear", restart_CS)) then
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        visc%Kv_shear(i,j,k) = Z2_T_rescale * visc%Kv_shear(i,j,k)
      enddo ; enddo ; enddo
    endif ; endif

    if (associated(visc%Kv_shear_Bu)) then ; if (query_initialized(visc%Kv_shear_Bu, "Kv_shear_Bu", restart_CS)) then
      do k=1,nz+1 ; do J=js-1,je ; do I=is-1,ie
        visc%Kv_shear_Bu(I,J,k) = Z2_T_rescale * visc%Kv_shear_Bu(I,J,k)
      enddo ; enddo ; enddo
    endif ; endif
  endif

  if (MLE_use_PBL_MLD .and. (Z_rescale /= 1.0)) then
    if (associated(visc%MLD)) then ; if (query_initialized(visc%MLD, "MLD", restart_CS)) then
      do j=js,je ; do i=is,ie
        visc%MLD(i,j) = Z_rescale * visc%MLD(i,j)
      enddo ; enddo
    endif ; endif
  endif

end subroutine set_visc_init

!> This subroutine dellocates any memory in the set_visc control structure.
subroutine set_visc_end(visc, CS)
  type(vertvisc_type), intent(inout) :: visc !< A structure containing vertical viscosities and
                                             !! related fields.  Elements are deallocated here.
  type(set_visc_CS),   pointer       :: CS   !< The control structure returned by a previous
                                             !! call to set_visc_init.
  if (CS%bottomdraglaw) then
    deallocate(visc%bbl_thick_u) ; deallocate(visc%bbl_thick_v)
    deallocate(visc%kv_bbl_u) ; deallocate(visc%kv_bbl_v)
    if (allocated(CS%bbl_u)) deallocate(CS%bbl_u)
    if (allocated(CS%bbl_v)) deallocate(CS%bbl_v)
  endif
  if (CS%Channel_drag) then
    deallocate(visc%Ray_u) ; deallocate(visc%Ray_v)
  endif
  if (CS%dynamic_viscous_ML) then
    deallocate(visc%nkml_visc_u) ; deallocate(visc%nkml_visc_v)
  endif
  if (associated(visc%Kd_shear)) deallocate(visc%Kd_shear)
  if (associated(visc%Kv_slow)) deallocate(visc%Kv_slow)
  if (associated(visc%TKE_turb)) deallocate(visc%TKE_turb)
  if (associated(visc%Kv_shear)) deallocate(visc%Kv_shear)
  if (associated(visc%Kv_shear_Bu)) deallocate(visc%Kv_shear_Bu)
  if (associated(visc%ustar_bbl)) deallocate(visc%ustar_bbl)
  if (associated(visc%TKE_bbl)) deallocate(visc%TKE_bbl)
  if (associated(visc%taux_shelf)) deallocate(visc%taux_shelf)
  if (associated(visc%tauy_shelf)) deallocate(visc%tauy_shelf)
  if (associated(visc%tbl_thick_shelf_u)) deallocate(visc%tbl_thick_shelf_u)
  if (associated(visc%tbl_thick_shelf_v)) deallocate(visc%tbl_thick_shelf_v)
  if (associated(visc%kv_tbl_shelf_u)) deallocate(visc%kv_tbl_shelf_u)
  if (associated(visc%kv_tbl_shelf_v)) deallocate(visc%kv_tbl_shelf_v)

  deallocate(CS)
end subroutine set_visc_end

!> \namespace mom_set_visc
!!
!! This would also be the module in which other viscous quantities that are flow-independent might be set.
!! This information is transmitted to other modules via a vertvisc type structure.
!!
!! The same code is used for the two velocity components, by indirectly referencing the velocities and
!! defining a handful of direction-specific defined variables.

end module MOM_set_visc
