!> Regrid columns for the SLight coordinate
module coord_slight

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_EOS,           only : EOS_type, calculate_compress
use MOM_EOS,           only : calculate_density, calculate_density_derivs
use regrid_interp,     only : interp_CS_type, regridding_set_ppolys
use regrid_interp,     only : NR_ITERATIONS, NR_TOLERANCE, DEGREE_MAX

implicit none ; private

!> Control structure containing required parameters for the SLight coordinate
type, public :: slight_CS ; private

  !> Number of layers/levels
  integer :: nk

  !> Minimum thickness allowed when building the new grid through regridding [H ~> m or kg m-2]
  real :: min_thickness

  !> Reference pressure for potential density calculations [R L2 T-2 ~> Pa]
  real :: ref_pressure

  !> Fraction (between 0 and 1) of compressibility to add to potential density
  !! profiles when interpolating for target grid positions. [nondim]
  real :: compressibility_fraction

  ! The following 4 parameters were introduced for use with the SLight coordinate:
  !> Depth over which to average to determine the mixed layer potential density [H ~> m or kg m-2]
  real :: Rho_ML_avg_depth

  !> Number of layers to offset the mixed layer density to find resolved stratification [nondim]
  real :: nlay_ml_offset

  !> The number of fixed-thickness layers at the top of the model
  integer :: nz_fixed_surface = 2

  !> The fixed resolution in the topmost SLight_nkml_min layers [H ~> m or kg m-2]
  real :: dz_ml_min

  !> If true, detect regions with much weaker stratification in the coordinate
  !! than based on in-situ density, and use a stretched coordinate there.
  logical :: fix_haloclines = .false.

  !> A length scale over which to filter T & S when looking for spuriously
  !! unstable water mass profiles [H ~> m or kg m-2].
  real :: halocline_filter_length

  !> A value of the stratification ratio that defines a problematic halocline region [nondim].
  real :: halocline_strat_tol

  !> Nominal density of interfaces [R ~> kg m-3].
  real, allocatable, dimension(:) :: target_density

  !> Maximum depths of interfaces [H ~> m or kg m-2].
  real, allocatable, dimension(:) :: max_interface_depths

  !> Maximum thicknesses of layers [H ~> m or kg m-2].
  real, allocatable, dimension(:) :: max_layer_thickness

  !> Interpolation control structure
  type(interp_CS_type) :: interp_CS
end type slight_CS

public init_coord_slight, set_slight_params, build_slight_column, end_coord_slight

contains

!> Initialise a slight_CS with pointers to parameters
subroutine init_coord_slight(CS, nk, ref_pressure, target_density, interp_CS, m_to_H)
  type(slight_CS),      pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,              intent(in) :: nk !< Number of layers in the grid
  real,                 intent(in) :: ref_pressure !< Coordinate reference pressure [R L2 T-2 ~> Pa]
  real, dimension(:),   intent(in) :: target_density !< Nominal density of interfaces [R ~> kg m-3]
  type(interp_CS_type), intent(in) :: interp_CS !< Controls for interpolation
  real,       optional, intent(in) :: m_to_H !< A conversion factor from m to the units of thicknesses

  real :: m_to_H_rescale  ! A unit conversion factor.

  if (associated(CS)) call MOM_error(FATAL, "init_coord_slight: CS already associated!")
  allocate(CS)
  allocate(CS%target_density(nk+1))

  m_to_H_rescale = 1.0 ; if (present(m_to_H)) m_to_H_rescale = m_to_H

  CS%nk                = nk
  CS%ref_pressure      = ref_pressure
  CS%target_density(:) = target_density(:)
  CS%interp_CS         = interp_CS

  ! Set real parameter default values
  CS%compressibility_fraction = 0. ! Nondim.
  CS%Rho_ML_avg_depth = 1.0 * m_to_H_rescale
  CS%nlay_ml_offset = 2.0          ! Nondim.
  CS%dz_ml_min = 1.0 * m_to_H_rescale
  CS%halocline_filter_length = 2.0 * m_to_H_rescale
  CS%halocline_strat_tol = 0.25    ! Nondim.

end subroutine init_coord_slight

!> This subroutine deallocates memory in the control structure for the coord_slight module
subroutine end_coord_slight(CS)
  type(slight_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
  if (.not. associated(CS)) return
  deallocate(CS%target_density)
  deallocate(CS)
end subroutine end_coord_slight

!> This subroutine can be used to set the parameters for the coord_slight module
subroutine set_slight_params(CS, max_interface_depths, max_layer_thickness, &
               min_thickness, compressibility_fraction, dz_ml_min, &
               nz_fixed_surface, Rho_ML_avg_depth, nlay_ML_offset, fix_haloclines, &
               halocline_filter_length, halocline_strat_tol, interp_CS)
  type(slight_CS),   pointer    :: CS !< Coordinate control structure
  real, dimension(:), &
           optional, intent(in) :: max_interface_depths !< Maximum depths of interfaces [H ~> m or kg m-2]
  real, dimension(:), &
           optional, intent(in) :: max_layer_thickness  !< Maximum thicknesses of layers [H ~> m or kg m-2]
  real,    optional, intent(in) :: min_thickness    !< Minimum thickness allowed when building the
                                      !! new grid through regridding [H ~> m or kg m-2]
  real,    optional, intent(in) :: compressibility_fraction !< Fraction (between 0 and 1) of
                                      !! compressibility to add to potential density profiles when
                                      !! interpolating for target grid positions. [nondim]
  real,    optional, intent(in) :: dz_ml_min        !< The fixed resolution in the topmost
                                      !! SLight_nkml_min layers [H ~> m or kg m-2]
  integer, optional, intent(in) :: nz_fixed_surface !< The number of fixed-thickness layers at the
                                      !! top of the model
  real,    optional, intent(in) :: Rho_ML_avg_depth !< Depth over which to average to determine
                                      !! the mixed layer potential density [H ~> m or kg m-2]
  real,    optional, intent(in) :: nlay_ML_offset   !< Number of layers to offset the mixed layer
                                      !! density to find resolved stratification [nondim]
  logical, optional, intent(in) :: fix_haloclines   !< If true, detect regions with much weaker than
                                      !! based on in-situ density, and use a stretched coordinate there.
  real,    optional, intent(in) :: halocline_filter_length !< A length scale over which to filter T & S
                                      !! when looking for spuriously unstable water mass profiles [H ~> m or kg m-2].
  real,    optional, intent(in) :: halocline_strat_tol !< A value of the stratification ratio that
                                      !! defines a problematic halocline region [nondim].
  type(interp_CS_type), &
           optional, intent(in) :: interp_CS !< Controls for interpolation

  if (.not. associated(CS)) call MOM_error(FATAL, "set_slight_params: CS not associated")

  if (present(max_interface_depths)) then
    if (size(max_interface_depths) /= CS%nk+1) &
      call MOM_error(FATAL, "set_slight_params: max_interface_depths inconsistent size")
    allocate(CS%max_interface_depths(CS%nk+1))
    CS%max_interface_depths(:) = max_interface_depths(:)
  endif

  if (present(max_layer_thickness)) then
    if (size(max_layer_thickness) /= CS%nk) &
      call MOM_error(FATAL, "set_slight_params: max_layer_thickness inconsistent size")
    allocate(CS%max_layer_thickness(CS%nk))
    CS%max_layer_thickness(:) = max_layer_thickness(:)
  endif

  if (present(min_thickness)) CS%min_thickness = min_thickness
  if (present(compressibility_fraction)) CS%compressibility_fraction = compressibility_fraction

  if (present(dz_ml_min)) CS%dz_ml_min = dz_ml_min
  if (present(nz_fixed_surface)) CS%nz_fixed_surface = nz_fixed_surface
  if (present(Rho_ML_avg_depth)) CS%Rho_ML_avg_depth = Rho_ML_avg_depth
  if (present(nlay_ML_offset)) CS%nlay_ML_offset = nlay_ML_offset
  if (present(fix_haloclines)) CS%fix_haloclines = fix_haloclines
  if (present(halocline_filter_length)) CS%halocline_filter_length = halocline_filter_length
  if (present(halocline_strat_tol)) then
    if (halocline_strat_tol > 1.0) call MOM_error(FATAL, "set_slight_params: "//&
        "HALOCLINE_STRAT_TOL must not exceed 1.0.")
    CS%halocline_strat_tol = halocline_strat_tol
  endif

  if (present(interp_CS)) CS%interp_CS = interp_CS
end subroutine set_slight_params

!> Build a SLight coordinate column
subroutine build_slight_column(CS, eqn_of_state, H_to_pres, H_subroundoff, &
                               nz, depth, h_col, T_col, S_col, p_col, z_col, z_col_new, &
                               h_neglect, h_neglect_edge)
  type(slight_CS),       intent(in)    :: CS    !< Coordinate control structure
  type(EOS_type),        pointer       :: eqn_of_state !< Equation of state structure
  real,                  intent(in)    :: H_to_pres !< A conversion factor from thicknesses to
                                                !! scaled pressure [R L2 T-2 H-1 ~> Pa m-1 or Pa m2 kg-1]
  real,                  intent(in)    :: H_subroundoff !< GV%H_subroundoff
  integer,               intent(in)    :: nz    !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive [H ~> m or kg m-2])
  real, dimension(nz),   intent(in)    :: T_col !< T for column
  real, dimension(nz),   intent(in)    :: S_col !< S for column
  real, dimension(nz),   intent(in)    :: h_col !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz),   intent(in)    :: p_col !< Layer center pressure [R L2 T-2 ~> Pa]
  real, dimension(nz+1), intent(in)    :: z_col !< Interface positions relative to the surface [H ~> m or kg m-2]
  real, dimension(nz+1), intent(inout) :: z_col_new !< Absolute positions of interfaces [H ~> m or kg m-2]
  real,        optional, intent(in)    :: h_neglect !< A negligibly small width for the purpose of
                                                !! cell reconstructions [H ~> m or kg m-2].
  real,        optional, intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose
                                                !! of edge value calculations [H ~> m or kg m-2].
  ! Local variables
  real, dimension(nz) :: rho_col        ! Layer densities [R ~> kg m-3]
  real, dimension(nz) :: T_f, S_f       ! Filtered layer temperature [degC] and salinity [ppt]
  logical, dimension(nz+1) :: reliable  ! If true, this interface is in a reliable position.
  real, dimension(nz+1) :: T_int, S_int ! Temperature [degC] and salinity [ppt] interpolated to interfaces.
  real, dimension(nz+1) :: rho_tmp      ! A temporary density [R ~> kg m-3]
  real, dimension(nz+1) :: drho_dp      ! The partial derivative of density with pressure [T2 L-2 ~> kg m-3 Pa-1]
  real, dimension(nz+1) :: p_IS, p_R    ! Pressures [R L2 T-2 ~> Pa]
  real, dimension(nz+1) :: drhoIS_dT    ! The partial derivative of in situ density with temperature
                                        ! in [R degC-1 ~> kg m-3 degC-1]
  real, dimension(nz+1) :: drhoIS_dS    ! The partial derivative of in situ density with salinity
                                        ! in [R ppt-1 ~> kg m-3 ppt-1]
  real, dimension(nz+1) :: drhoR_dT     ! The partial derivative of reference density with temperature
                                        ! in [R degC-1 ~> kg m-3 degC-1]
  real, dimension(nz+1) :: drhoR_dS     ! The partial derivative of reference density with salinity
                                        ! in [R ppt-1 ~> kg m-3 ppt-1]
  real, dimension(nz+1) :: strat_rat
  real :: H_to_cPa    ! A conversion factor from thicknesses to the compressibility fraction times
                      ! the units of pressure [R L2 T-2 H-1 ~> Pa m-1 or Pa m2 kg-1]
  real :: drIS, drR   ! In situ and reference density differences [R ~> kg m-3]
  real :: Fn_now, I_HStol, Fn_zero_val ! Nondimensional variables [nondim]
  real :: z_int_unst  ! The depth where the stratification allows the interior grid to start [H ~> m or kg m-2]
  real :: dz          ! A uniform layer thickness in very shallow water [H ~> m or kg m-2].
  real :: dz_ur       ! The total thickness of an unstable region [H ~> m or kg m-2].
  real :: wgt, cowgt  ! A weight and its complement [nondim].
  real :: rho_ml_av   ! The average potential density in a near-surface region [R ~> kg m-3].
  real :: H_ml_av     ! A thickness to try to use in taking the near-surface average [H ~> m or kg m-2].
  real :: rho_x_z     ! A cumulative integral of a density [R H ~> kg m-2 or kg2 m-5].
  real :: z_wt        ! The thickness actually used in taking the near-surface average [H ~> m or kg m-2].
  real :: k_interior  ! The (real) value of k where the interior grid starts [nondim].
  real :: k_int2      ! The (real) value of k where the interior grid starts [nondim].
  real :: z_interior  ! The depth where the interior grid starts [H ~> m or kg m-2].
  real :: z_ml_fix    ! The depth at which the fixed-thickness near-surface layers end [H ~> m or kg m-2].
  real :: dz_dk       ! The thickness of layers between the fixed-thickness
                      ! near-surface layars and the interior [H ~> m or kg m-2].
  real :: Lfilt       ! A filtering lengthscale [H ~> m or kg m-2].
  logical :: maximum_depths_set ! If true, the maximum depths of interface have been set.
  logical :: maximum_h_set      ! If true, the maximum layer thicknesses have been set.
  real :: k2_used, k2here, dz_sum, z_max
  integer :: k2
  real :: h_tr, b_denom_1, b1, d1 ! Temporary variables used by the tridiagonal solver.
  real, dimension(nz) :: c1  ! Temporary variables used by the tridiagonal solver.
  integer :: kur1, kur2  ! The indicies at the top and bottom of an unreliable region.
  integer :: kur_ss      ! The index to start with in the search for the next unstable region.
  integer :: i, j, k, nkml

  maximum_depths_set = allocated(CS%max_interface_depths)
  maximum_h_set = allocated(CS%max_layer_thickness)

  if (z_col(nz+1) - z_col(1) < nz*CS%min_thickness) then
    ! This is a nearly massless total depth, so distribute the water evenly.
    dz = (z_col(nz+1) - z_col(1)) / real(nz)
    do K=2,nz ; z_col_new(K) = z_col(1) + dz*real(K-1) ; enddo
  else
    call calculate_density(T_col, S_col, p_col, rho_col, eqn_of_state)

    ! Find the locations of the target potential densities, flagging
    ! locations in apparently unstable regions as not reliable.
    call rho_interfaces_col(rho_col, h_col, z_col, CS%target_density, nz, &
                            z_col_new, CS, reliable, debug=.true., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)

    ! Ensure that the interfaces are at least CS%min_thickness apart.
    if (CS%min_thickness > 0.0) then
      ! Move down interfaces below overly thin layers.
      do K=2,nz ; if (z_col_new(K) < z_col_new(K-1) + CS%min_thickness) then
        z_col_new(K) = z_col_new(K-1) + CS%min_thickness
      endif ; enddo
      ! Now move up any interfaces that are too close to the bottom.
      do K=nz,2,-1 ; if (z_col_new(K) > z_col_new(K+1) - CS%min_thickness) then
        z_col_new(K) = z_col_new(K+1) - CS%min_thickness
      else
        exit ! No more interfaces can be too close to the bottom.
      endif ; enddo
    endif

    ! Fix up the unreliable regions.
    kur_ss = 2 ! reliable(1) and reliable(nz+1) must always be true.
    do
      ! Search for the uppermost unreliable interface postion.
      kur1 = nz+2
      do K=kur_ss,nz ; if (.not.reliable(K)) then
        kur1 = K ; exit
      endif ; enddo
      if (kur1 > nz) exit ! Everything is now reliable.

      kur2 = kur1-1 ! For error checking.
      do K=kur1+1,nz+1 ; if (reliable(K)) then
        kur2 = K-1 ; kur_ss = K ; exit
      endif ; enddo
      if (kur2 < kur1) call MOM_error(FATAL, "Bad unreliable range.")

      dz_ur = z_col_new(kur2+1) - z_col_new(kur1-1)
  !        drho = CS%target_density(kur2+1) - CS%target_density(kur1-1)
      ! Perhaps reset the wgt and cowgt depending on how bad the old interface
      ! locations were.
      wgt = 1.0 ; cowgt = 0.0 ! = 1.0-wgt
      do K=kur1,kur2
        z_col_new(K) = cowgt*z_col_new(K) + &
              wgt * (z_col_new(kur1-1) + dz_ur*(K - (kur1-1)) / ((kur2 - kur1) + 2))
      enddo
    enddo

    ! Determine which interfaces are in the s-space region and the depth extent
    ! of this region.
    z_wt = 0.0 ; rho_x_z = 0.0
    H_ml_av = CS%Rho_ml_avg_depth
    do k=1,nz
      if (z_wt + h_col(k) >= H_ml_av) then
        rho_x_z = rho_x_z + rho_col(k) * (H_ml_av - z_wt)
        z_wt = H_ml_av
        exit
      else
        rho_x_z =  rho_x_z + rho_col(k) * h_col(k)
        z_wt = z_wt + h_col(k)
      endif
    enddo
    if (z_wt > 0.0) rho_ml_av = rho_x_z / z_wt

    nkml = CS%nz_fixed_surface
    ! Find the interface that matches rho_ml_av.
    if (rho_ml_av <= CS%target_density(nkml)) then
      k_interior = CS%nlay_ml_offset + real(nkml)
    elseif (rho_ml_av > CS%target_density(nz+1)) then
      k_interior = real(nz+1)
    else ; do K=nkml,nz
      if ((rho_ml_av >= CS%target_density(K)) .and. &
          (rho_ml_av <  CS%target_density(K+1))) then
        k_interior = (CS%nlay_ml_offset + K) + &
                (rho_ml_av - CS%target_density(K)) / &
                (CS%target_density(K+1) - CS%target_density(K))
        exit
      endif
    enddo ; endif
    if (k_interior > real(nz+1)) k_interior = real(nz+1)

    ! Linearly interpolate to find z_interior.  This could be made more sophisticated.
    K = int(ceiling(k_interior))
    z_interior = (K-k_interior)*z_col_new(K-1) + (1.0+(k_interior-K))*z_col_new(K)

    if (CS%fix_haloclines) then
  !       ! Identify regions above the reference pressure where the chosen
  !       ! potential density significantly underestimates the actual
  !       ! stratification, and use these to find a second estimate of
  !       ! z_int_unst and k_interior.

      if (CS%halocline_filter_length > 0.0) then
        Lfilt = CS%halocline_filter_length

        ! Filter the temperature and salnity with a fixed lengthscale.
        h_tr = h_col(1) + H_subroundoff
        b1 = 1.0 / (h_tr + Lfilt) ; d1 = h_tr * b1
        T_f(1) = (b1*h_tr)*T_col(1) ;  S_f(1) = (b1*h_tr)*S_col(1)
        do k=2,nz
          c1(k) = Lfilt * b1
          h_tr = h_col(k) + H_subroundoff ; b_denom_1 = h_tr + d1*Lfilt
          b1 = 1.0 / (b_denom_1 + Lfilt) ; d1 = b_denom_1 * b1
          T_f(k) = b1 * (h_tr*T_col(k) + Lfilt*T_f(k-1))
          S_f(k) = b1 * (h_tr*S_col(k) + Lfilt*S_f(k-1))
        enddo
        do k=nz-1,1,-1
          T_f(k) = T_f(k) + c1(k+1)*T_f(k+1) ; S_f(k) = S_f(k) + c1(k+1)*S_f(k+1)
        enddo
      else
        do k=1,nz ; T_f(k) = T_col(k) ; S_f(k) = S_col(k) ; enddo
      endif

      T_int(1) = T_f(1) ; S_int(1) = S_f(1)
      do K=2,nz
        T_int(K) = 0.5*(T_f(k-1) + T_f(k)) ; S_int(K) = 0.5*(S_f(k-1) + S_f(k))
        p_IS(K) = z_col(K) * H_to_pres
        p_R(K) = CS%ref_pressure + CS%compressibility_fraction * ( p_IS(K) - CS%ref_pressure )
      enddo
      T_int(nz+1) = T_f(nz) ; S_int(nz+1) = S_f(nz)
      p_IS(nz+1) = z_col(nz+1) * H_to_pres
      call calculate_density_derivs(T_int, S_int, p_IS, drhoIS_dT, drhoIS_dS, &
                                    eqn_of_state, (/2,nz/) )
      call calculate_density_derivs(T_int, S_int, p_R, drhoR_dT, drhoR_dS, &
                                    eqn_of_state, (/2,nz/) )
      if (CS%compressibility_fraction > 0.0) then
        call calculate_compress(T_int, S_int, p_R(:), rho_tmp, drho_dp, 2, nz-1, eqn_of_state)
      else
        do K=2,nz ; drho_dp(K) = 0.0 ; enddo
      endif

      H_to_cPa = CS%compressibility_fraction * H_to_pres
      strat_rat(1) = 1.0
      do K=2,nz
        drIS = drhoIS_dT(K) * (T_f(k) - T_f(k-1)) + &
               drhoIS_dS(K) * (S_f(k) - S_f(k-1))
        drR = (drhoR_dT(K) * (T_f(k) - T_f(k-1)) + &
               drhoR_dS(K) * (S_f(k) - S_f(k-1))) + &
              drho_dp(K) * (H_to_cPa*0.5*(h_col(k) + h_col(k-1)))

        if (drIS <= 0.0) then
          strat_rat(K) = 2.0 ! Maybe do this? => ; if (drR < 0.0) strat_rat(K) = -2.0
        else
          strat_rat(K) = 2.0*max(drR,0.0) / (drIS + abs(drR))
        endif
      enddo
      strat_rat(nz+1) = 1.0

      z_int_unst = 0.0 ; Fn_now = 0.0
      Fn_zero_val = min(2.0*CS%halocline_strat_tol, &
                        0.5*(1.0 + CS%halocline_strat_tol))
      if (CS%halocline_strat_tol > 0.0) then
        ! Use Adcroft's reciprocal rule.
        I_HStol = 0.0 ; if (Fn_zero_val - CS%halocline_strat_tol > 0.0) &
          I_HStol = 1.0 / (Fn_zero_val - CS%halocline_strat_tol)
        do k=nz,1,-1 ; if (CS%ref_pressure > p_IS(k+1)) then
          z_int_unst = z_int_unst + Fn_now * h_col(k)
          if (strat_rat(K) <= Fn_zero_val) then
            if (strat_rat(K) <= CS%halocline_strat_tol) then ; Fn_now = 1.0
            else
              Fn_now = max(Fn_now, (Fn_zero_val - strat_rat(K)) * I_HStol)
            endif
          endif
        endif ; enddo
      else
        do k=nz,1,-1 ; if (CS%ref_pressure > p_IS(k+1)) then
          z_int_unst = z_int_unst + Fn_now * h_col(k)
          if (strat_rat(K) <= CS%halocline_strat_tol) Fn_now = 1.0
        endif ; enddo
      endif

      if (z_interior < z_int_unst) then
        ! Find a second estimate of the extent of the s-coordinate region.
        kur1 = max(int(ceiling(k_interior)),2)
        if (z_col_new(kur1-1) < z_interior) then
          k_int2 = kur1
          do K = kur1,nz+1 ; if (z_col_new(K) >= z_int_unst) then
            ! This is linear interpolation again.
            if (z_col_new(K-1) >= z_int_unst) &
              call MOM_error(FATAL,"build_grid_SLight, bad halocline structure.")
            k_int2 = real(K-1) + (z_int_unst - z_col_new(K-1)) / &
                                     (z_col_new(K) - z_col_new(K-1))
            exit
          endif ; enddo
          if (z_col_new(nz+1) < z_int_unst) then
            ! This should be unnecessary.
            z_int_unst = z_col_new(nz+1) ; k_int2 = real(nz+1)
          endif

          ! Now take the larger values.
          if (k_int2 > k_interior) then
            k_interior = k_int2 ; z_interior = z_int_unst
          endif
        endif
      endif
    endif  ! fix_haloclines

    z_col_new(1) = 0.0
    do K=2,nkml+1
      z_col_new(K) = min((K-1)*CS%dz_ml_min, &
                         z_col_new(nz+1) - CS%min_thickness*(nz+1-K))
    enddo
    z_ml_fix = z_col_new(nkml+1)
    if (z_interior > z_ml_fix) then
      dz_dk = (z_interior - z_ml_fix) / (k_interior - (nkml+1))
      do K=nkml+2,int(floor(k_interior))
        z_col_new(K) = z_ml_fix + dz_dk * (K - (nkml+1))
      enddo
    else ! The fixed-thickness z-region penetrates into the interior.
      do K=nkml+2,nz
        if (z_col_new(K) <= z_col_new(CS%nz_fixed_surface+1)) then
          z_col_new(K) = z_col_new(CS%nz_fixed_surface+1)
        else ; exit ; endif
      enddo
    endif

    if (maximum_depths_set .and. maximum_h_set) then ; do k=2,nz
      ! The loop bounds are 2 & nz so the top and bottom interfaces do not move.
      ! Recall that z_col_new is positive downward.
      z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K), &
                         z_col_new(K-1) + CS%max_layer_thickness(k-1))
    enddo ; elseif (maximum_depths_set) then ; do K=2,nz
      z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K))
    enddo ; elseif (maximum_h_set) then ; do k=2,nz
      z_col_new(K) = min(z_col_new(K), z_col_new(K-1) + CS%max_layer_thickness(k-1))
    enddo ; endif

  endif ! Total thickness exceeds nz*CS%min_thickness.

end subroutine build_slight_column

!> Finds the new interface locations in a column of water that match the
!! prescribed target densities.
subroutine rho_interfaces_col(rho_col, h_col, z_col, rho_tgt, nz, z_col_new, &
                              CS, reliable, debug, h_neglect, h_neglect_edge)
  integer,               intent(in)    :: nz      !< Number of layers
  real, dimension(nz),   intent(in)    :: rho_col !< Initial layer reference densities [R ~> kg m-3].
  real, dimension(nz),   intent(in)    :: h_col   !< Initial layer thicknesses [H ~> m or kg m-2].
  real, dimension(nz+1), intent(in)    :: z_col   !< Initial interface heights [H ~> m or kg m-2].
  real, dimension(nz+1), intent(in)    :: rho_tgt !< Interface target densities.
  real, dimension(nz+1), intent(inout) :: z_col_new !< New interface heights [H ~> m or kg m-2].
  type(slight_CS),       intent(in)    :: CS      !< Coordinate control structure
  logical, dimension(nz+1), intent(inout) :: reliable !< If true, the interface positions
                                                  !! are well defined from a stable region.
  logical,     optional, intent(in)    :: debug   !< If present and true, do debugging checks.
  real,        optional, intent(in)    :: h_neglect !< A negligibly small width for the purpose of
                                                  !! cell reconstructions [H ~> m or kg m-2]
  real,        optional, intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose
                                                  !! of edge value calculations [H ~> m or kg m-2]

  real, dimension(nz+1) :: ru_max_int ! The maximum and minimum densities in
  real, dimension(nz+1) :: ru_min_int ! an unstable region around an interface [R ~> kg m-3].
  real, dimension(nz)   :: ru_max_lay ! The maximum and minimum densities in
  real, dimension(nz)   :: ru_min_lay ! an unstable region containing a layer [R ~> kg m-3].
  real, dimension(nz,2) :: ppoly_i_E  ! Edge value of polynomial [R ~> kg m-3]
  real, dimension(nz,2) :: ppoly_i_S  ! Edge slope of polynomial [R H-1 ~> kg m-4 or m-1]
  real, dimension(nz,DEGREE_MAX+1) :: ppoly_i_coefficients ! Coefficients of polynomial [R ~> kg m-3]
  logical, dimension(nz)   :: unstable_lay ! If true, this layer is in an unstable region.
  logical, dimension(nz+1) :: unstable_int ! If true, this interface is in an unstable region.
  real :: rt  ! The current target density [R ~> kg m-3].
  real :: zf  ! The fractional z-position within a layer of the target density [nondim].
  real :: rfn ! The target density relative to the interpolated density [R ~> kg m-3]
  real :: a(5) ! Coefficients of a local polynomial minus the target density [R ~> kg m-3].
  real :: zf1, zf2   ! Two previous estimates of zf [nondim]
  real :: rfn1, rfn2 ! Values of rfn at zf1 and zf2 [R ~> kg m-3]
  real :: drfn_dzf   ! The partial derivative of rfn with zf [R ~> kg m-3]
  real :: sgn, delta_zf, zf_prev ! [nondim]
  real :: tol  ! The tolerance for convergence of zf [nondim]
  logical :: k_found ! If true, the position has been found.
  integer :: k_layer ! The index of the stable layer containing an interface.
  integer :: ppoly_degree
  integer :: k, k1, k1_min, itt, max_itt, m

  real :: z_sgn  ! 1 or -1, depending on whether z increases with increasing K.
  logical :: debugging

  debugging = .false. ; if (present(debug)) debugging = debug
  max_itt = NR_ITERATIONS
  tol = NR_TOLERANCE

  z_sgn = 1.0 ; if ( z_col(1) > z_col(nz+1) ) z_sgn = -1.0
  if (debugging) then
    do K=1,nz
      if (abs((z_col(K+1) - z_col(K)) - z_sgn*h_col(k)) > &
          1.0e-14*(abs(z_col(K+1)) + abs(z_col(K)) + abs(h_col(k))) ) &
        call MOM_error(FATAL, "rho_interfaces_col: Inconsistent z_col and h_col")
    enddo
  endif

  if ( z_col(1) == z_col(nz+1) ) then
    ! This is a massless column!
    do K=1,nz+1 ; z_col_new(K) = z_col(1) ; reliable(K) = .true. ; enddo
    return
  endif

  ! This sets up the piecewise polynomials based on the rho_col profile.
  call regridding_set_ppolys(CS%interp_CS, rho_col, nz, h_col, ppoly_i_E, ppoly_i_S, &
       ppoly_i_coefficients, ppoly_degree, h_neglect, h_neglect_edge)

  ! Determine the density ranges of unstably stratified segments.
  ! Interfaces that start out in an unstably stratified segment can
  ! only escape if they are outside of the bounds of that segment, and no
  ! interfaces are ever mapped into an unstable segment.
  unstable_int(1) = .false.
  ru_max_int(1) = ppoly_i_E(1,1)

  unstable_lay(1) = (ppoly_i_E(1,1) > ppoly_i_E(1,2))
  ru_max_lay(1) = max(ppoly_i_E(1,1), ppoly_i_E(1,2))

  do K=2,nz
    unstable_int(K) = (ppoly_i_E(k-1,2) > ppoly_i_E(k,1))
    ru_max_int(K) = max(ppoly_i_E(k-1,2), ppoly_i_E(k,1))
    ru_min_int(K) = min(ppoly_i_E(k-1,2), ppoly_i_E(k,1))
    if (unstable_int(K) .and. unstable_lay(k-1)) &
      ru_max_int(K) = max(ru_max_lay(k-1), ru_max_int(K))

    unstable_lay(k) = (ppoly_i_E(k,1) > ppoly_i_E(k,2))
    ru_max_lay(k) = max(ppoly_i_E(k,1), ppoly_i_E(k,2))
    ru_min_lay(k) = min(ppoly_i_E(k,1), ppoly_i_E(k,2))
    if (unstable_lay(k) .and. unstable_int(K)) &
      ru_max_lay(k) = max(ru_max_int(K), ru_max_lay(k))
  enddo
  unstable_int(nz+1) = .false.
  ru_min_int(nz+1) = ppoly_i_E(nz,2)

  do K=nz,1,-1
    if (unstable_lay(k) .and. unstable_int(K+1)) &
      ru_min_lay(k) = min(ru_min_int(K+1), ru_min_lay(k))

    if (unstable_int(K) .and. unstable_lay(k)) &
      ru_min_int(K) = min(ru_min_lay(k), ru_min_int(K))
  enddo

  z_col_new(1) = z_col(1) ; reliable(1) = .true.
  k1_min = 1
  do K=2,nz ! Find the locations of the various target densities for the interfaces.
    rt = rho_tgt(K)
    k_layer = -1
    k_found = .false.

    ! Many light layers are found at the top, so start there.
    if (rt <= ppoly_i_E(k1_min,1)) then
      z_col_new(K) = z_col(k1_min)
      k_found = .true.
      ! Do not change k1_min for the next layer.
    elseif (k1_min == nz+1) then
      z_col_new(K) = z_col(nz+1)
    else
      ! Start with the previous location and search outward.
      if (unstable_int(K) .and. (rt >= ru_min_int(K)) .and. (rt <= ru_max_int(K))) then
        ! This interface started in an unstable region and should not move due to remapping.
        z_col_new(K) = z_col(K) ; reliable(K) = .false.
        k1_min = K ; k_found = .true.
      elseif ((rt >= ppoly_i_E(k-1,2)) .and. (rt <= ppoly_i_E(k,1))) then
        ! This interface is already in the right place and does not move.
        z_col_new(K) = z_col(K) ; reliable(K) = .true.
        k1_min = K ; k_found = .true.
      elseif (rt < ppoly_i_E(k-1,2)) then   ! Search upward
        do k1=K-1,k1_min,-1
          ! Check whether rt is in layer k.
          if ((rt < ppoly_i_E(k1,2)) .and. (rt > ppoly_i_E(k1,1))) then
            ! rt is in layer k.
            k_layer = k1
            k1_min = k1 ; k_found = .true. ; exit
          elseif (unstable_lay(k1) .and. (rt >= ru_min_lay(k1)) .and. (rt <= ru_max_lay(K1))) then
            ! rt would be found at unstable layer that it can not penetrate.
            !   It is possible that this can never happen?
            z_col_new(K) = z_col(K1+1) ; reliable(K) = .false.
            k1_min = k1 ; k_found = .true. ; exit
          endif
          ! Check whether rt is at interface K.
          if (k1 > 1) then ; if ((rt <= ppoly_i_E(k1,1)) .and. (rt >= ppoly_i_E(k1-1,2))) then
            ! rt is at interface K1
            z_col_new(K) = z_col(K1) ; reliable(K) = .true.
            k1_min = k1 ; k_found = .true. ; exit
          elseif (unstable_int(K1) .and. (rt >= ru_min_int(k1)) .and. (rt <= ru_max_int(K1))) then
            ! rt would be found at an unstable interface that it can not pass.
            !   It is possible that this can never happen?
            z_col_new(K) = z_col(K1) ; reliable(K) = .false.
            k1_min = k1 ; k_found = .true. ; exit
          endif ; endif
        enddo

        if (.not.k_found) then
          ! This should not happen unless k1_min = 1.
          if (k1_min < 2) then
            z_col_new(K) = z_col(k1_min)
          else
            z_col_new(K) = z_col(k1_min)
          endif
        endif

      else  ! Search downward
        do k1=K,nz
          if ((rt < ppoly_i_E(k1,2)) .and. (rt > ppoly_i_E(k1,1))) then
            ! rt is in layer k.
            k_layer = k1
            k1_min = k1 ; k_found = .true. ; exit
          elseif (unstable_lay(k1) .and. (rt >= ru_min_lay(k1)) .and. (rt <= ru_max_lay(K1))) then
            ! rt would be found at unstable layer that it can not penetrate.
            !   It is possible that this can never happen?
            z_col_new(K) = z_col(K1)
            reliable(K) = .false.
            k1_min = k1 ; k_found = .true. ; exit
          endif
          if (k1 < nz) then ; if ((rt <= ppoly_i_E(k1+1,1)) .and. (rt >= ppoly_i_E(k1,2))) then
            ! rt is at interface K1+1

            z_col_new(K) = z_col(K1+1) ; reliable(K) = .true.
            k1_min = k1+1 ; k_found = .true. ; exit
          elseif (unstable_int(K1+1) .and. (rt >= ru_min_int(k1+1)) .and. (rt <= ru_max_int(K1+1))) then
            ! rt would be found at an unstable interface that it can not pass.
            !   It is possible that this can never happen?
            z_col_new(K) = z_col(K1+1)
            reliable(K) = .false.
            k1_min = k1+1 ; k_found = .true. ; exit
          endif ; endif
        enddo
        if (.not.k_found) then
          z_col_new(K) = z_col(nz+1)
          if (rt >= ppoly_i_E(nz,2)) then
            reliable(K) = .true.
          else
            reliable(K) = .false.
          endif
        endif
      endif

      if (k_layer > 0) then  ! The new location is inside of layer k_layer.
        ! Note that this is coded assuming that this layer is stably stratified.
        if (.not.(ppoly_i_E(k1,2) > ppoly_i_E(k1,1))) call MOM_error(FATAL, &
          "build_grid_SLight: Erroneously searching for an interface in an unstratified layer.")

        ! Use the false position method to find the location (degree <= 1) or the first guess.
        zf = (rt - ppoly_i_E(k1,1)) / (ppoly_i_E(k1,2) - ppoly_i_E(k1,1))

        if (ppoly_degree > 1) then ! Iterate to find the solution.
          a(:) = 0.0 ; a(1) = ppoly_i_coefficients(k_layer,1) - rt
          do m=2,ppoly_degree+1 ; a(m) = ppoly_i_coefficients(k_layer,m) ; enddo
          ! Bracket the root.
          zf1 = 0.0 ; rfn1 = a(1)
          zf2 = 1.0 ; rfn2 =  a(1) + (a(2) + (a(3) + (a(4) + a(5))))
          if (rfn1 * rfn2 > 0.0) call MOM_error(FATAL, "build_grid_SLight: Bad bracketing.")

          do itt=1,max_itt
            rfn = a(1) + zf*(a(2) + zf*(a(3) + zf*(a(4) + zf*a(5))))
            ! Reset one of the ends of the bracket.
            if (rfn * rfn1 > 0.0) then
              zf1 = zf ; rfn1 = rfn
            else
              zf2 = zf ; rfn2 = rfn
            endif
            if (rfn1 == rfn2) exit

            drfn_dzf = (a(2) + zf*(2.0*a(3) + zf*(3.0*a(4) + zf*4.0*a(5))))
            sgn = 1.0 ; if (drfn_dzf < 0.0) sgn = -1.0

            if ((sgn*(zf - rfn) >= zf1 * abs(drfn_dzf)) .and. &
                (sgn*(zf - rfn) <= zf2 * abs(drfn_dzf))) then
              delta_zf = -rfn / drfn_dzf
              zf = zf + delta_zf
            else ! Newton's method goes out of bounds, so use a false position method estimate
              zf_prev = zf
              zf = ( rfn2 * zf1 - rfn1 * zf2 ) / (rfn2 - rfn1)
              delta_zf = zf - zf_prev
            endif

            if (abs(delta_zf) < tol) exit
          enddo
        endif
        z_col_new(K) = z_col(k_layer) + zf * z_sgn * h_col(k_layer)
        reliable(K) = .true.
      endif

    endif

  enddo
  z_col_new(nz+1) = z_col(nz+1) ; reliable(nz+1) = .true.

end subroutine rho_interfaces_col

end module coord_slight
