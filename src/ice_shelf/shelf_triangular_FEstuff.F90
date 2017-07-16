module shelf_triangular_FEstuff

use MOM_diag_mediator, only : diag_ctrl, time_type, enable_averaging, disable_averaging
use MOM_grid, only :  ocean_grid_type
use MOM_time_manager, only : time_type, set_time, time_type_to_real
use MOM_restart, only : restart_init, restore_state, MOM_restart_CS
use MOM_EOS, only : EOS_type
use user_shelf_init, only : user_ice_shelf_CS

implicit none ; private

#include <MOM_memory.h>
type, public :: ice_shelf_CS ; private
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()
  type(ocean_grid_type) :: grid ! A structure containing metrics, etc.
  ! The rest is private
  character(len=128) :: restart_output_dir = ' '
  real, pointer, dimension(:,:) :: &
    mass_shelf => NULL(), &   ! The mass per unit area of the ice shelf or sheet, in kg m-2.
    area_shelf_h => NULL(), & ! The area per cell covered by the ice shelf, in m2.

    t_flux => NULL(), &   ! The UPWARD sensible ocean heat flux at the ocean-ice
                    ! interface, in W m-2.
    salt_flux => NULL(), & ! The downward salt flux at the ocean-ice interface, in kg m-2 s-1.
    lprec => NULL(), &    ! The downward liquid water flux at the ocean-ice interface,
                    ! in kg m-2 s-1.
    ! Perhaps these diagnostics should only be kept with the call?
    exch_vel_t => NULL(), &
    exch_vel_s => NULL(), &
    tfreeze => NULL(), &  ! The freezing point potential temperature an the ice-ocean
                    ! interface, in deg C.
    tflux_shelf => NULL(), & ! The UPWARD diffusive heat flux in the ice shelf at the
                    ! ice-ocean interface, in W m-2.
!!! DNG !!!
    u_shelf => NULL(), & ! the zonal (?) velocity of the ice shelf/sheet... in meters per second???
          ! on q-points (B grid)
    v_shelf => NULL(), & ! the meridional velocity of the ice shelf/sheet... m/s ??
          ! on q-points (B grid)
    h_shelf => NULL(), & ! the thickness of the shelf in m... redundant with mass
          ! but may make code more readable
    hmask => NULL(),&    ! used to indicate ice-covered cells, as well as partially-covered
          ! 1: fully covered, solve for velocity here
          !   (for now all ice-covered cells are treated the same, this may change)
          ! 2: partially covered, do not solve for velocity
          ! 0: no ice in cell.
          ! 3: bdry condition on thickness set - not in computational domain
          ! -2 : default (out of computational boundary, and not = 3

          ! NOTE: hmask will change over time and NEEDS TO BE MAINTAINED
          !    otherwise the wrong nodes will be included in velocity calcs.
    u_face_mask => NULL(), v_face_mask => NULL(), &
   ! masks for velocity boundary conditions - on *C GRID* - this is because the FEM solution
   !      cares about FACES THAT GET INTEGRATED OVER, not vertices
   ! Will represent boundary conditions on computational boundary (or permanent boundary
   !      between fast-moving and near-stagnant ice
   ! FOR NOW: 1=interior bdry, 0=no-flow boundary, 2=stress bdry condition, 3=inhomogeneous dirichlet boundary
    umask => NULL(), vmask => NULL(), &
   ! masks on the actual degrees of freedom (B grid) -
   !   1=normal node, 3=inhomogeneous boundary node, 0 - no flow node (will also get ice-free nodes)
    ice_visc_bilinear => NULL(), &
    ice_visc_lower_tri => NULL(), &
    ice_visc_upper_tri => NULL(), &
    thickness_boundary_values => NULL(), &
    u_boundary_values => NULL(), &
    v_boundary_values => NULL(), &


    taub_beta_eff_bilinear => NULL(), & ! nonlinear part of "linearized" basal stress - exact form depends on basal law exponent
                ! and/or whether flow is "hybridized" a la Goldberg 2011
    taub_beta_eff_lower_tri => NULL(), &
    taub_beta_eff_upper_tri => NULL(), &

    OD_rt => NULL(), float_frac_rt => NULL(), &
    OD_av => NULL(), float_frac => NULL()  !! two arrays that represent averages of ocean values that are maintained
                  !! within the ice shelf module and updated based on the "ocean state".
                  !! OD_av is ocean depth, and float_frac is the average amount of time
                  !! a cell is "exposed", i.e. the column thickness is below a threshold.
                  !! both are averaged over the time of a diagnostic (ice velocity)

                       !! [if float_frac = 1 ==> grounded; obv. counterintuitive; might fix]

  real :: ustar_bg     ! A minimum value for ustar under ice shelves, in m s-1.
  real :: Cp           ! The heat capacity of sea water, in J kg-1 K-1.
  real :: Cp_ice       ! The heat capacity of fresh ice, in J kg-1 K-1.
  real :: gamma_t      !   The (fixed) turbulent exchange velocity in the
                       ! 2-equation formulation, in m s-1.
  real :: Salin_ice    ! The salinity of shelf ice, in PSU.
  real :: Temp_ice     ! The core temperature of shelf ice, in C.
  real :: kv_ice       ! The viscosity of ice, in m2 s-1.
  real :: density_ice  ! A typical density of ice, in kg m-3.
  real :: kv_molec     ! The molecular kinematic viscosity of sea water, m2 s-1.
  real :: kd_molec_salt  ! The molecular diffusivity of salt, in m2 s-1.
  real :: kd_molec_temp  ! The molecular diffusivity of heat, in m2 s-1.
  real :: Lat_fusion   ! The latent heat of fusion, in J kg-1.

!!!! PHYSICAL AND NUMERICAL PARAMETERS FOR ICE DYNAMICS !!!!!!

  real :: time_step    ! this is the shortest timestep that the ice shelf sees, and
            ! is equal to the forcing timestep (it is passed in when the shelf
            ! is initialized - so need to reorganize MOM driver.
            ! it will be the prognistic timestep ... maybe.

!!! all need to be initialized

  real :: A_glen_isothermal
  real :: n_glen
  real :: eps_glen_min
  real :: C_basal_friction
  real :: n_basal_friction
  real :: density_ocean_avg    ! this does not affect ocean circulation OR thermodynamics
                ! it is to estimate the gravitational driving force at the shelf front
                ! (until we think of a better way to do it- but any difference will be negligible)
  real :: thresh_float_col_depth ! the water column depth over which the shelf if considered to be floating
  real :: input_flux
  real :: input_thickness

  real :: len_lat ! this really should be a Grid or Domain field


  real :: velocity_update_time_step ! the time to update the velocity through the nonlinear
                    ! elliptic equation. i think this should be done no more often than
                    ! ~ once a day (maybe longer) because it will depend on ocean values
                    ! that are averaged over this time interval, and the solve will begin
                    ! to lose meaning if it is done too frequently
  integer :: velocity_update_sub_counter ! there is no outer loop for the velocity solve; the counter will have to be stored
  integer :: velocity_update_counter ! the "outer" timestep number
  integer :: nstep_velocity        ! ~ (velocity_update_time_step / time_step)

  real :: cg_tolerance, nonlinear_tolerance
  integer :: cg_max_iterations
  integer :: nonlin_solve_err_mode  ! 1: exit vel solve based on nonlin residual
                    ! 2: exit based on "fixed point" metric (|u - u_last| / |u| < tol where | | is infty-norm
  real    :: CFL_factor            ! in uncoupled run, how to limit subcycled advective timestep
                      ! i.e. dt = CFL_factor * min (dx / u)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(time_type) :: Time ! The component's time.
  type(EOS_type), pointer :: eqn_of_state => NULL() ! Type that indicates the
                                        ! equation of state to use.
  logical :: isshelf   ! True if a shelf model is to be used.
  logical :: shelf_mass_is_dynamic ! True if the ice shelf mass changes with
                       ! time.
  logical :: override_shelf_movement ! If true, user code specifies the shelf
                       ! movement instead of using the dynamic ice-shelf mode.
  logical :: isthermo  ! True if the ice shelf can exchange heat and mass with
                       ! the underlying ocean.
  logical :: threeeq   ! If true, the 3 equation consistency equations are
                       ! used to calculate the flux at the ocean-ice interface.
  integer :: id_melt = -1, id_exch_vel_s = -1, id_exch_vel_t = -1, &
             id_tfreeze = -1, id_tfl_shelf = -1, &
         id_u_shelf = -1, id_v_shelf = -1, id_h_shelf = -1, id_h_mask = -1, &
         id_u_mask = -1, id_v_mask = -1, &
         id_surf_elev = -1, id_bathym = -1, id_float_frac = -1, id_col_thick = -1, &
         id_area_shelf_h = -1, id_OD_rt = -1, id_float_frac_rt = -1
  type(diag_ctrl) :: diag     ! A structure that is used to control diagnostic
                              ! output.
  type(user_ice_shelf_CS), pointer :: user_CS => NULL()

  logical :: write_output_to_file ! this is for seeing arrays w/out netcdf capability
end type ice_shelf_CS
contains

subroutine matrix_diagonal_triangle (CS, u_diagonal, v_diagonal)

  type(ice_shelf_CS),    pointer       :: CS
  real, dimension (:,:), intent(inout) :: u_diagonal, v_diagonal

! returns the diagonal entries of the matrix for a Jacobi preconditioning

  real, pointer, dimension (:,:)       :: umask, vmask, &
                          nu_lower, nu_upper, beta_lower, beta_upper, hmask
  type(ocean_grid_type), pointer :: G
  integer ::  i, j, is, js, cnt, isc, jsc, iec, jec
  real :: A, n, ux, uy, vx, vy, eps_min, domain_width, dxh, dyh, dxdyh

  G => CS%grid

!   if (G%symmetric) then
!     isym=1
!   else
!     isym=0
!   endif


  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  umask => CS%umask ; vmask => CS%vmask ; hmask => CS%hmask
  nu_lower => CS%ice_visc_lower_tri ; nu_upper => CS%ice_visc_upper_tri
  beta_lower => CS%taub_beta_eff_lower_tri ; beta_upper => CS%taub_beta_eff_upper_tri

  do i=isc-1,iec+1  ; do j=jsc-1,jec+1 ; if (hmask(i,j) .eq. 1) then
    dxh = G%dxT(i,j)
    dyh = G%dyT(i,j)
    dxdyh = G%areaT(i,j)

    if (umask (i,j-1) .eq. 1) then ! this (bot right) is a degree of freedom node

      ux = 1./dxh ; uy = 0./dyh
      vx = 0. ; vy = 0.

      u_diagonal (i,j-1) = u_diagonal (i,j-1) + &
          .5 * dxdyh * nu_lower (i,j) * ((4*ux+2*vy) * (1./dxh) + (uy+vy) * (0./dyh))

      u_diagonal (i,j-1) = u_diagonal (i,j-1) + &
          beta_lower(i,j) * dxdyh * 1./24

      ux = 0. ; uy = 0.
      vx = 1./dxh ; vy = 0./dyh

      v_diagonal (i,j-1) = v_diagonal (i,j-1) + &
          .5 * dxdyh * nu_lower (i,j) * ((uy+vx) * (1./dxh) + (4*vy+2*ux) * (0./dyh))

      v_diagonal (i,j-1) = v_diagonal (i,j-1) + &
          beta_lower(i,j) * dxdyh * 1./24

      ux = 0./dxh ; uy = -1./dyh
      vx = 0. ; vy = 0.

      u_diagonal (i,j-1) = u_diagonal (i,j-1) + &
          .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (0./dxh) + (uy+vy) * (-1./dyh))

      u_diagonal (i,j-1) = u_diagonal (i,j-1) + &
          beta_upper(i,j) * dxdyh * 1./24

      vx = 0./dxh ; vy = -1./dyh
      ux = 0. ; uy = 0.

      v_diagonal (i,j-1) = v_diagonal (i,j-1) + &
          .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (0./dxh) + (4*vy+2*ux) * (-1./dyh))

      v_diagonal (i,j-1) = v_diagonal (i,j-1) + &
          beta_upper(i,j) * dxdyh * 1./24

    endif

    if (umask (i-1,j) .eq. 1) then ! this (top left) is a degree of freedom node

      ux = 0./dxh ; uy = 1./dyh
      vx = 0. ; vy = 0.

      u_diagonal (i-1,j) = u_diagonal (i-1,j) + &
          .5 * dxdyh * nu_lower (i,j) * ((4*ux+2*vy) * (0./dxh) + (uy+vy) * (1./dyh))

      u_diagonal (i,j-1) = u_diagonal (i,j-1) + &
          beta_lower(i,j) * dxdyh * 1./24

      ux = 0. ; uy = 0.
      vx = 0./dxh ; vy = 1./dyh

      v_diagonal (i-1,j) = v_diagonal (i-1,j) + &
          .5 * dxdyh * nu_lower (i,j) * ((uy+vx) * (0./dxh) + (4*vy+2*ux) * (1./dyh))

      v_diagonal (i,j-1) = v_diagonal (i,j-1) + &
          beta_lower(i,j) * dxdyh * 1./24

      ux = -1./dxh ; uy = 0./dyh
      vx = 0. ; vy = 0.

      u_diagonal (i-1,j) = u_diagonal (i-1,j) + &
          .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (-1./dxh) + (uy+vy) * (0./dyh))

      u_diagonal (i,j-1) = u_diagonal (i,j-1) + &
          beta_upper(i,j) * dxdyh * 1./24

      vx = -1./dxh ; vy = 0./dyh
      ux = 0. ; uy = 0.

      v_diagonal (i-1,j) = v_diagonal (i-1,j) + &
          .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (-1./dxh) + (4*vy+2*ux) * (0./dyh))

      v_diagonal (i,j-1) = v_diagonal (i,j-1) + &
          beta_upper(i,j) * dxdyh * 1./24

    endif

    if (umask (i-1,j-1) .eq. 1) then ! this (bot left) is a degree of freedom node

      ux = -1./dxh ; uy = -1./dyh
      vx = 0. ; vy = 0.

      u_diagonal (i-1,j-1) = u_diagonal (i-1,j-1) + &
          .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (-1./dxh) + (uy+vy) * (-1./dyh))

      u_diagonal (i-1,j-1) = u_diagonal (i-1,j-1) + &
          beta_lower(i,j) * dxdyh * 1./24

      vx = -1./dxh ; vy = -1./dyh
      ux = 0. ; uy = 0.

      v_diagonal (i-1,j-1) = v_diagonal (i-1,j-1) + &
          .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (-1./dxh) + (4*vy+2*ux) * (-1./dyh))

      v_diagonal (i-1,j-1) = v_diagonal (i-1,j-1) + &
          beta_lower(i,j) * dxdyh * 1./24
    endif

    if (umask (i,j) .eq. 1) then ! this (top right) is a degree of freedom node

      ux = 1./ dxh ; uy = 1./dyh
      vx = 0. ; vy = 0.

      u_diagonal (i,j) = u_diagonal (i,j) + &
          .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (1./dxh) + (uy+vy) * (1./dyh))

      u_diagonal (i,j) = u_diagonal (i,j) + &
          beta_upper(i,j) * dxdyh * 1./24

      vx = 1./ dxh ; vy = 1./dyh
      ux = 0. ; uy = 0.

      v_diagonal (i,j) = v_diagonal (i,j) + &
          .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (1./dxh) + (4*vy+2*ux) * (1./dyh))

      v_diagonal (i,j) = v_diagonal (i,j) + &
          beta_upper(i,j) * dxdyh * 1./24

    endif
  endif ; enddo ; enddo

end subroutine matrix_diagonal_triangle

!~ subroutine apply_boundary_values_triangle (CS, time, u_boundary_contr, v_boundary_contr)

  !~ type(time_type),       intent(in)    :: Time
  !~ type(ice_shelf_CS),    pointer       :: CS
  !~ real, dimension (:,:), intent(inout) :: u_boundary_contr, v_boundary_contr

!~ ! this will be a per-setup function. the boundary values of thickness and velocity
!~ ! (and possibly other variables) will be updated in this function

  !~ real, pointer, dimension (:,:)       :: u_boundary_values, &
                          !~ v_boundary_values, &
                          !~ umask, vmask, hmask, &
                          !~ nu_lower, nu_upper, beta_lower, beta_upper
  !~ type(ocean_grid_type), pointer :: G
  !~ integer :: 0, i, j, cnt, isc, jsc, iec, jec
  !~ real :: A, n, ux, uy, vx, vy, eps_min, domain_width, dxh, dyh, dxdyh

  !~ G => CS%grid

!~ !   if (G%symmetric) then
!~ !    isym=1
!~ !   else
!~ !    isym=0
!~ !   endif



  !~ isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  !~ u_boundary_values => CS%u_boundary_values
  !~ v_boundary_values => CS%v_boundary_values
  !~ umask => CS%umask ; vmask => CS%vmask ; hmask => CS%hmask
  !~ nu_lower => CS%ice_visc_lower_tri ; nu_upper => CS%ice_visc_upper_tri
  !~ beta_lower => CS%taub_beta_eff_lower_tri ; beta_upper => CS%taub_beta_eff_upper_tri

  !~ domain_width = CS%len_lat

  !~ do i=isc-1,iec+1 ; do j=jsc-1,jec+1 ; if (hmask(i,j) .eq. 1) then

    !~ if ((umask(i-1,j-1) .eq. 3) .OR. (umask(i,j-1) .eq. 3) .OR. (umask(i-1,j) .eq. 3)) then

      !~ dxh = G%dxh(i,j)
      !~ dyh = G%dyh(i,j)
      !~ dxdyh = G%dxdyh(i,j)

      !~ ux = (u_boundary_values(i,j-1)-u_boundary_values(i-1,j-1))/dxh
      !~ vx = (v_boundary_values(i,j-1)-v_boundary_values(i-1,j-1))/dxh
      !~ uy = (u_boundary_values(i-1,j)-u_boundary_values(i-1,j-1))/dyh
      !~ vy = (v_boundary_values(i-1,j)-v_boundary_values(i-1,j-1))/dyh

      !~ if (umask (i,j-1) .eq. 1) then ! this (bot right) is a degree of freedom node

        !~ u_boundary_contr (i,j-1) = u_boundary_contr (i,j-1) + &
            !~ .5 * dxdyh * nu_lower (i,j) * ((4*ux+2*vy) * (1./dxh) + (uy+vy) * (0./dyh))

        !~ v_boundary_contr (i,j-1) = v_boundary_contr (i,j-1) + &
            !~ .5 * dxdyh * nu_lower (i,j) * ((uy+vx) * (1./dxh) + (4*vy+2*ux) * (0./dyh))

        !~ u_boundary_contr (i,j-1) = u_boundary_contr (i,j-1) + &
            !~ beta_lower(i,j) * dxdyh * 1./24 * (u_boundary_values(i-1,j-1) + &
                          !~ u_boundary_values(i-1,j) + u_boundary_values(i,j-1))

        !~ v_boundary_contr (i,j-1) = v_boundary_contr (i,j-1) + &
            !~ beta_lower(i,j) * dxdyh * 1./24 * (v_boundary_values(i-1,j-1) + &
                         !~ v_boundary_values(i-1,j) + v_boundary_values(i,j-1))
      !~ endif

      !~ if (umask (i-1,j) .eq. 1) then ! this (top left) is a degree of freedom node

        !~ u_boundary_contr (i-1,j) = u_boundary_contr (i-1,j) + &
            !~ .5 * dxdyh * nu_lower (i,j) * ((4*ux+2*vy) * (0./dxh) + (uy+vy) * (1./dyh))

        !~ v_boundary_contr (i-1,j) = v_boundary_contr (i-1,j) + &
            !~ .5 * dxdyh * nu_lower (i,j) * ((uy+vx) * (0./dxh) + (4*vy+2*ux) * (1./dyh))

        !~ u_boundary_contr (i,j-1) = u_boundary_contr (i,j-1) + &
            !~ beta_lower(i,j) * dxdyh * 1./24 * (u_boundary_values(i-1,j-1) + &
                            !~ u_boundary_values(i-1,j) + u_boundary_values(i,j-1))

        !~ v_boundary_contr (i,j-1) = v_boundary_contr (i,j-1) + &
            !~ beta_lower(i,j) * dxdyh * 1./24 * (v_boundary_values(i-1,j-1) + &
                            !~ v_boundary_values(i-1,j) + v_boundary_values(i,j-1))
      !~ endif

      !~ if (umask (i-1,j-1) .eq. 1) then ! this (bot left) is a degree of freedom node

        !~ u_boundary_contr (i-1,j-1) = u_boundary_contr (i-1,j-1) + &
            !~ .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (-1./dxh) + (uy+vy) * (-1./dyh))

        !~ v_boundary_contr (i-1,j-1) = v_boundary_contr (i-1,j-1) + &
            !~ .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (-1./dxh) + (4*vy+2*ux) * (-1./dyh))

        !~ u_boundary_contr (i-1,j-1) = u_boundary_contr (i-1,j-1) + &
            !~ beta_lower(i,j) * dxdyh * 1./24 * (u_boundary_values(i-1,j-1) + &
                            !~ u_boundary_values(i-1,j) + u_boundary_values(i,j-1))

        !~ v_boundary_contr (i-1,j-1) = v_boundary_contr (i-1,j-1) + &
            !~ beta_lower(i,j) * dxdyh * 1./24 * (v_boundary_values(i-1,j-1) + &
                            !~ v_boundary_values(i-1,j) + v_boundary_values(i,j-1))
      !~ endif

    !~ endif

    !~ if ((umask(i,j) .eq. 3) .OR. (umask(i,j-1) .eq. 3) .OR. (umask(i-1,j) .eq. 3)) then

      !~ dxh = G%dxh(i,j)
      !~ dyh = G%dyh(i,j)
      !~ dxdyh = G%dxdyh(i,j)

      !~ ux = (u_boundary_values(i,j)-u_boundary_values(i-1,j))/dxh
      !~ vx = (v_boundary_values(i,j)-v_boundary_values(i-1,j))/dxh
      !~ uy = (u_boundary_values(i,j)-u_boundary_values(i,j-1))/dyh
      !~ vy = (v_boundary_values(i,j)-v_boundary_values(i,j-1))/dyh

      !~ if (umask (i,j-1) .eq. 1) then ! this (bot right) is a degree of freedom node

          !~ u_boundary_contr (i,j-1) = u_boundary_contr (i,j-1) + &
              !~ .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (0./dxh) + (uy+vy) * (-1./dyh))

          !~ v_boundary_contr (i,j-1) = v_boundary_contr (i,j-1) + &
              !~ .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (0./dxh) + (4*vy+2*ux) * (-1./dyh))

          !~ u_boundary_contr (i,j-1) = u_boundary_contr (i,j-1) + &
              !~ beta_upper(i,j) * dxdyh * 1./24 * (u_boundary_values(i,j) + &
                                      !~ u_boundary_values(i-1,j) +  &
                                !~ u_boundary_values(i,j-1))

          !~ v_boundary_contr (i,j-1) = v_boundary_contr (i,j-1) + &
              !~ beta_upper(i,j) * dxdyh * 1./24 * (u_boundary_values(i,j) + &
                                      !~ u_boundary_values(i-1,j) +  &
                                !~ u_boundary_values(i,j-1))
      !~ endif

      !~ if (umask (i-1,j) .eq. 1) then ! this (top left) is a degree of freedom node

        !~ u_boundary_contr (i-1,j) = u_boundary_contr (i-1,j) + &
            !~ .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (-1./dxh) + (uy+vy) * (0./dyh))

        !~ v_boundary_contr (i-1,j) = v_boundary_contr (i-1,j) + &
            !~ .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (-1./dxh) + (4*vy+2*ux) * (0./dyh))

        !~ u_boundary_contr (i,j-1) = u_boundary_contr (i,j-1) + &
            !~ beta_upper(i,j) * dxdyh * 1./24 * (u_boundary_values(i,j) + &
                                    !~ u_boundary_values(i-1,j) +  &
                              !~ u_boundary_values(i,j-1))

        !~ v_boundary_contr (i,j-1) = v_boundary_contr (i,j-1) + &
            !~ beta_upper(i,j) * dxdyh * 1./24 * (u_boundary_values(i,j) + &
                                    !~ u_boundary_values(i-1,j) +  &
                              !~ u_boundary_values(i,j-1))
      !~ endif

      !~ if (umask (i,j) .eq. 1) then ! this (top right) is a degree of freedom node

        !~ u_boundary_contr (i,j) = u_boundary_contr (i,j) + &
            !~ .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (1./dxh) + (uy+vy) * (1./dyh))

        !~ v_boundary_contr (i,j) = v_boundary_contr (i,j) + &
            !~ .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (1./dxh) + (4*vy+2*ux) * (1./dyh))

        !~ u_boundary_contr (i,j) = u_boundary_contr (i,j) + &
            !~ beta_upper(i,j) * dxdyh * 1./24 * (u_boundary_values(i,j) + &
                                    !~ u_boundary_values(i-1,j) +  &
                              !~ u_boundary_values(i,j-1))

        !~ v_boundary_contr (i,j) = v_boundary_contr (i,j) + &
            !~ beta_upper(i,j) * dxdyh * 1./24 * (u_boundary_values(i,j) + &
                                    !~ u_boundary_values(i-1,j) +  &
                              !~ u_boundary_values(i,j-1))
      !~ endif


    !~ endif
  !~ endif ; enddo ; enddo

!~ end subroutine apply_boundary_values_triangle

!~ subroutine calc_shelf_visc_triangular (CS,u,v)
  !~ type(ice_shelf_CS),         pointer   :: CS
  !~ real, dimension(:,:), intent(inout)    :: u, v

!~ ! update DEPTH_INTEGRATED viscosity, based on horizontal strain rates - this is for triangle FEM solve so there is
!~ ! an "upper" and "lower" triangular viscosity

!~ ! also this subroutine updates the nonlinear part of the basal traction

!~ ! this may be subject to change later... to make it "hybrid"

  !~ real, pointer, dimension (:,:)    :: nu_lower , &
                         !~ nu_upper, &
                       !~ beta_eff_lower, &
                       !~ beta_eff_upper
  !~ real, pointer, dimension (:,:)    :: H,    &! thickness
                       !~ hmask

  !~ type(ocean_grid_type), pointer :: G
  !~ integer :: 0, i, j, iscq, iecq, jscq, jecq, isd, jsd, ied, jed, iegq, jegq, giec, gjec, gisc, gjsc, cnt, isc, jsc, iec, jec, is, js
  !~ real :: A, n, ux, uy, vx, vy, eps_min, umid, vmid, unorm, C_basal_friction, n_basal_friction, dxh, dyh, dxdyh

  !~ G => CS%grid

  !~ isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  !~ iscq = G%iscq ; iecq = G%iecq ; jscq = G%jscq ; jecq = G%jecq
  !~ isd = G%isd ; jsd = G%jsd ; ied = G%isd ; jed = G%jsd
  !~ iegq = G%iegq ; jegq = G%jegq
  !~ gisc = G%domain%nx_halo+1 ; gjsc = G%domain%ny_halo+1
  !~ giec = G%domain%nxtot+gisc ; gjec = G%domain%nytot+gjsc
  !~ is = iscq - (1-0); js = jscq - (1-0)

  !~ A = CS%A_glen_isothermal ; n = CS%n_glen; eps_min = CS%eps_glen_min

  !~ H => CS%h_shelf
  !~ hmask => CS%hmask
  !~ nu_upper => CS%ice_visc_upper_tri
  !~ nu_lower => CS%ice_visc_lower_tri
  !~ beta_eff_upper => CS%taub_beta_eff_upper_tri
  !~ beta_eff_lower => CS%taub_beta_eff_lower_tri

  !~ C_basal_friction = CS%C_basal_friction ; n_basal_friction = CS%n_basal_friction

  !~ do i=isd,ied
    !~ do j=jsd,jed

      !~ dxh = G%dxh(i,j)
      !~ dyh = G%dyh(i,j)
      !~ dxdyh = G%dxdyh(i,j)

      !~ if (hmask (i,j) .eq. 1) then
        !~ ux = (u(i,j-1)-u(i-1,j-1)) / dxh
        !~ vx = (v(i,j-1)-v(i-1,j-1)) / dxh
        !~ uy = (u(i-1,j)-u(i-1,j-1)) / dyh
        !~ vy = (v(i-1,j)-v(i-1,j-1)) / dyh

        !~ nu_lower(i,j) = A**(-1/n) * (ux**2+vy**2+ux*vy.25*(uy+vx)**2+eps_min**2) ** ((1-n)/(2*n)) * H(i,j)
        !~ umid = 1./3 * (u(i-1,j-1)+u(i-1,j)+u(i,j-1))
        !~ vmid = 1./3 * (v(i-1,j-1)+v(i-1,j)+v(i,j-1))
        !~ unorm = sqrt (umid**2+vmid**2+(eps_min*dxh)**2) ; beta_eff_lower (i,j) = C_basal_friction * unorm ** (n_basal_friction-1)

        !~ ux = (u(i,j)-u(i-1,j)) / dxh
        !~ vx = (v(i,j)-v(i-1,j)) / dxh
        !~ uy = (u(i,j)-u(i,j-1)) / dyh
        !~ vy = (u(i,j)-u(i,j-1)) / dyh

        !~ nu_upper(i,j) = A**(-1/n) * (ux**2+vy**2+ux*vy.25*(uy+vx)**2+eps_min**2) ** ((1-n)/(2*n)) * H(i,j)
        !~ umid = 1./3 * (u(i,j)+u(i-1,j)+u(i,j-1))
        !~ vmid = 1./3 * (v(i,j)+v(i-1,j)+v(i,j-1))
        !~ unorm = sqrt (umid**2+vmid**2+(eps_min*dxh)**2) ; beta_eff_upper (i,j) = C_basal_friction * unorm ** (n_basal_friction-1)

      !~ endif
    !~ enddo
  !~ enddo

!~ end subroutine calc_shelf_visc_triangular


!~ subroutine CG_action_triangular (uret, vret, u, v, umask, vmask, hmask, nu_upper, nu_lower, &
                !~ beta_upper, beta_lower, dxh, dyh, dxdyh, is, ie, js, je, 0)

!~ real, dimension (:,:), intent (inout)  :: uret, vret
!~ real, dimension (:,:), intent (in)     :: u, v
!~ real, dimension (:,:), intent (in)     :: umask, vmask
!~ real, dimension (:,:), intent (in)     :: hmask, nu_upper, nu_lower, beta_upper, beta_lower
!~ real, dimension (:,:), intent (in)     :: dxh, dyh, dxdyh
!~ integer, intent(in)               :: is, ie, js, je, 0

!~ ! the linear action of the matrix on (u,v) with triangular finite elements
!~ ! as of now everything is passed in so no grid pointers or anything of the sort have to be dereferenced,
!~ ! but this may change pursuant to conversations with others
!~ !
!~ ! is & ie are the cells over which the iteration is done; this may change between calls to this subroutine
!~ !     in order to make less frequent halo updates
!~ ! isym = 1 if grid is symmetric, 0 o.w.

  !~ real :: ux, uy, vx, vy
  !~ integer :: i,j

  !~ do i=is,ie
    !~ do j=js,je

      !~ if (hmask(i,j) .eq. 1) then ! this cell's vertices contain degrees of freedom

        !~ ux = (u(i,j-1)-u(i-1,j-1))/dxh(i,j)
        !~ vx = (v(i,j-1)-v(i-1,j-1))/dxh(i,j)
        !~ uy = (u(i-1,j)-u(i-1,j-1))/dyh(i,j)
        !~ vy = (v(i-1,j)-v(i-1,j-1))/dyh(i,j)

        !~ if (umask(i,j-1) .eq. 1) then ! this (bot right) is a degree of freedom node

          !~ uret(i,j-1) = uret(i,j-1) + &
              !~ .5 * dxdyh(i,j) * nu_lower (i,j) * ((4*ux+2*vy) * (1./dxh(i,j)) + (uy+vy) * (0./dyh(i,j)))

          !~ vret(i,j-1) = vret(i,j-1) + &
              !~ .5 * dxdyh(i,j) * nu_lower (i,j) * ((uy+vx) * (1./dxh(i,j)) + (4*vy+2*ux) * (0./dyh(i,j)))

          !~ uret(i,j-1) = uret(i,j-1) + &
              !~ beta_lower(i,j) * dxdyh(i,j) * 1./24 * (u(i-1,j-1) + &
                                      !~ u(i-1,j) + u(i,j-1))

          !~ vret(i,j-1) = vret(i,j-1) + &
              !~ beta_lower(i,j) * dxdyh(i,j) * 1./24 * (v(i-1,j-1) + &
                                      !~ v(i-1,j) + v(i,j-1))
        !~ endif

        !~ if (umask(i-1,j) .eq. 1) then ! this (top left) is a degree of freedom node

          !~ uret(i-1,j) = uret(i-1,j) + &
              !~ .5 * dxdyh(i,j) * nu_lower (i,j) * ((4*ux+2*vy) * (0./dxh(i,j)) + (uy+vy) * (1./dyh(i,j)))

          !~ vret(i-1,j) = vret(i-1,j) + &
              !~ .5 * dxdyh(i,j) * nu_lower (i,j) * ((uy+vx) * (0./dxh(i,j)) + (4*vy+2*ux) * (1./dyh(i,j)))

          !~ uret(i,j-1) = uret(i,j-1) + &
              !~ beta_lower(i,j) * dxdyh(i,j) * 1./24 * (u(i-1,j-1) + &
                                      !~ u(i-1,j) + u(i,j-1))

          !~ vret(i,j-1) = vret(i,j-1) + &
              !~ beta_lower(i,j) * dxdyh(i,j) * 1./24 * (v(i-1,j-1) + &
                                      !~ v(i-1,j) + v(i,j-1))
        !~ endif

        !~ if (umask(i-1,j-1) .eq. 1) then ! this (bot left) is a degree of freedom node

          !~ uret(i-1,j-1) = uret(i-1,j-1) + &
              !~ .5 * dxdyh(i,j) * nu_upper (i,j) * ((4*ux+2*vy) * (-1./dxh(i,j)) + (uy+vy) * (-1./dyh(i,j)))

          !~ vret(i-1,j-1) = vret(i-1,j-1) + &
              !~ .5 * dxdyh(i,j) * nu_upper (i,j) * ((uy+vx) * (-1./dxh(i,j)) + (4*vy+2*ux) * (-1./dyh(i,j)))

          !~ uret(i-1,j-1) = uret(i-1,j-1) + &
              !~ beta_lower(i,j) * dxdyh(i,j) * 1./24 * (u(i-1,j-1) + &
                                      !~ u(i-1,j) + u(i,j-1))

          !~ vret(i-1,j-1) = vret(i-1,j-1) + &
              !~ beta_lower(i,j) * dxdyh(i,j) * 1./24 * (v(i-1,j-1) + &
                                      !~ v(i-1,j) + v(i,j-1))
        !~ endif


        !~ ux = (u(i,j)-u(i-1,j))/dxh(i,j)
        !~ vx = (v(i,j)-v(i-1,j))/dxh(i,j)
        !~ uy = (u(i,j)-u(i,j-1))/dyh(i,j)
        !~ vy = (v(i,j)-v(i,j-1))/dyh(i,j)

        !~ if (umask(i,j-1) .eq. 1) then ! this (bot right) is a degree of freedom node

          !~ uret(i,j-1) = uret(i,j-1) + &
              !~ .5 * dxdyh(i,j) * nu_upper (i,j) * ((4*ux+2*vy) * (0./dxh(i,j)) + (uy+vy) * (-1./dyh(i,j)))

          !~ vret(i,j-1) = vret(i,j-1) + &
              !~ .5 * dxdyh(i,j) * nu_upper (i,j) * ((uy+vx) * (0./dxh(i,j)) + (4*vy+2*ux) * (-1./dyh(i,j)))

          !~ uret(i,j-1) = uret(i,j-1) + &
              !~ beta_upper(i,j) * dxdyh(i,j) * 1./24 * (u(i,j) + &
                                      !~ u(i-1,j) + u(i,j-1))

          !~ vret(i,j-1) = vret(i,j-1) + &
              !~ beta_upper(i,j) * dxdyh(i,j) * 1./24 * (u(i,j) + &
                                      !~ u(i-1,j) + u(i,j-1))
        !~ endif

        !~ if (umask(i-1,j) .eq. 1) then ! this (top left) is a degree of freedom node

          !~ uret(i-1,j) = uret(i-1,j) + &
              !~ .5 * dxdyh(i,j) * nu_upper (i,j) * ((4*ux+2*vy) * (-1./dxh(i,j)) + (uy+vy) * (0./dyh(i,j)))

          !~ vret(i-1,j) = vret(i-1,j) + &
              !~ .5 * dxdyh(i,j) * nu_upper (i,j) * ((uy+vx) * (-1./dxh(i,j)) + (4*vy+2*ux) * (0./dyh(i,j)))

          !~ uret(i,j-1) = uret(i,j-1) + &
              !~ beta_upper(i,j) * dxdyh(i,j) * 1./24 * (u(i,j) + &
                                      !~ u(i-1,j) + u(i,j-1))

          !~ vret(i,j-1) = vret(i,j-1) + &
              !~ beta_upper(i,j) * dxdyh(i,j) * 1./24 * (u(i,j) + &
                                      !~ u(i-1,j) + u(i,j-1))
        !~ endif

        !~ if (umask(i,j) .eq. 1) then ! this (top right) is a degree of freedom node

          !~ uret(i,j) = uret(i,j) + &
              !~ .5 * dxdyh(i,j) * nu_upper (i,j) * ((4*ux+2*vy) * (1./dxh(i,j)) + (uy+vy) * (1./dyh(i,j)))

          !~ vret(i,j) = vret(i,j) + &
              !~ .5 * dxdyh(i,j) * nu_upper (i,j) * ((uy+vx) * (1./dxh(i,j)) + (4*vy+2*ux) * (1./dyh(i,j)))

          !~ uret(i,j) = uret(i,j) + &
              !~ beta_upper(i,j) * dxdyh(i,j) * 1./24 * (u(i,j) + &
                                      !~ u(i-1,j) + u(i,j-1))

          !~ vret(i,j) = vret(i,j) + &
              !~ beta_upper(i,j) * dxdyh(i,j) * 1./24 * (u(i,j) + &
                                      !~ u(i-1,j) + u(i,j-1))
        !~ endif

      !~ endif

    !~ enddo
  !~ enddo

!~ end subroutine CG_action_triangular


END MODULE shelf_triangular_FEstuff
