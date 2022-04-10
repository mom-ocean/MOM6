!> This module contains the hybgen unmixing routines from HYCOM, with
!! modifications to follow the MOM6 coding conventions and several bugs fixed
module MOM_hybgen_unmix

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_EOS,             only : EOS_type, calculate_density, calculate_density_derivs
use MOM_error_handler,   only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, param_file_type, log_param
use MOM_hybgen_regrid,   only : hybgen_column_init
use MOM_hybgen_regrid,   only : hybgen_regrid_CS, get_hybgen_regrid_params
use MOM_tracer_registry, only : tracer_registry_type, tracer_type, MOM_tracer_chkinv
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!> Control structure containing required parameters for the hybgen coordinate generator
type, public :: hybgen_unmix_CS ; private

  integer :: nsigma  !< Number of sigma levels used by HYBGEN
  real :: hybiso     !< Hybgen uses PCM if layer is within hybiso of target density [R ~> kg m-3]

  real :: dp00i   !< Deep isopycnal spacing minimum thickness [H ~> m or kg m-2]
  real :: qhybrlx !< Hybgen relaxation amount per thermodynamic time steps [nondim]

  real, allocatable, dimension(:) ::  &
    dp0k, &     !< minimum deep    z-layer separation [H ~> m or kg m-2]
    ds0k        !< minimum shallow z-layer separation [H ~> m or kg m-2]

  real :: dpns  !< depth to start terrain following [H ~> m or kg m-2]
  real :: dsns  !< depth to stop terrain following [H ~> m or kg m-2]
  real :: min_dilate !< The minimum amount of dilation that is permitted when converting target
                     !! coordinates from z to z* [nondim].  This limit applies when wetting occurs.
  real :: max_dilate !< The maximum amount of dilation that is permitted when converting target
                     !! coordinates from z to z* [nondim].  This limit applies when drying occurs.

  real :: topiso_const !< Shallowest depth for isopycnal layers [H ~> m or kg m-2]
  ! real, dimension(:,:), allocatable :: topiso

  real :: ref_pressure !< Reference pressure for density calculations [R L2 T-2 ~> Pa]
  real, allocatable, dimension(:) :: target_density !< Nominal density of interfaces [R ~> kg m-3]

end type hybgen_unmix_CS

public hybgen_unmix, init_hybgen_unmix, end_hybgen_unmix
public set_hybgen_unmix_params

contains

!> Initialise a hybgen_unmix_CS control structure and store its parameters
subroutine init_hybgen_unmix(CS, GV, US, param_file, hybgen_regridCS)
  type(hybgen_unmix_CS),   pointer    :: CS  !< Unassociated pointer to hold the control structure
  type(verticalGrid_type), intent(in) :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file
  type(hybgen_regrid_CS),  pointer    :: hybgen_regridCS !< Control structure for hybgen
                                             !! regridding for sharing parameters.

  character(len=40)               :: mdl = "MOM_hybgen" ! This module's name.
  integer :: k

  if (associated(CS)) call MOM_error(FATAL, "init_hybgen_unmix: CS already associated!")
  allocate(CS)
  allocate(CS%target_density(GV%ke))

  allocate(CS%dp0k(GV%ke), source=0.0) ! minimum deep z-layer separation
  allocate(CS%ds0k(GV%ke), source=0.0) ! minimum shallow z-layer separation

  ! Set the parameters for the hybgen unmixing from a hybgen regridding control structure.
  call get_hybgen_regrid_params(hybgen_regridCS, ref_pressure=CS%ref_pressure, &
                nsigma=CS%nsigma, dp0k=CS%dp0k, ds0k=CS%ds0k, &
                dp00i=CS%dp00i, topiso_const=CS%topiso_const, qhybrlx=CS%qhybrlx, &
                hybiso=CS%hybiso, min_dilate=CS%min_dilate, max_dilate=CS%max_dilate, &
                target_density=CS%target_density)

  ! Determine the depth range over which to use a sigma (terrain-following) coordinate.
  ! --- terrain following starts at depth dpns and ends at depth dsns
  if (CS%nsigma == 0) then
    CS%dpns = CS%dp0k(1)
    CS%dsns = 0.0
  else
    CS%dpns = 0.0
    CS%dsns = 0.0
    do k=1,CS%nsigma
      CS%dpns = CS%dpns + CS%dp0k(k)
      CS%dsns = CS%dsns + CS%ds0k(k)
    enddo !k
  endif !nsigma

end subroutine init_hybgen_unmix

!> This subroutine deallocates memory in the control structure for the hybgen unmixing module
subroutine end_hybgen_unmix(CS)
  type(hybgen_unmix_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
  if (.not. associated(CS)) return

  deallocate(CS%target_density)
  deallocate(CS%dp0k, CS%ds0k)
  deallocate(CS)
end subroutine end_hybgen_unmix

!> This subroutine can be used to set the parameters for the hybgen module
subroutine set_hybgen_unmix_params(CS, min_thickness)
  type(hybgen_unmix_CS),  pointer    :: CS !< Coordinate unmixing control structure
  real,    optional, intent(in) :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]

  if (.not. associated(CS)) call MOM_error(FATAL, "set_hybgen_params: CS not associated")

!  if (present(min_thickness)) CS%min_thickness = min_thickness
end subroutine set_hybgen_unmix_params


!> Unmix the properties in the lowest layer with mass if it is too light, and make
!! any other changes to the water column to prepare for regridding.
subroutine hybgen_unmix(G, GV, US, CS, tv, Reg, ntr, h)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(hybgen_unmix_CS),   intent(in)    :: CS  !< hybgen control structure
  type(thermo_var_ptrs),   intent(inout) :: tv  !< Thermodynamics structure
  type(tracer_registry_type), pointer    :: Reg !< Tracer registry structure
  integer,                 intent(in)    :: ntr !< The number of tracers in the registry, or
                                                !! 0 if the registry is not in use.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h   !< Layer thicknesses [H ~> m or kg m-2]

! --- --------------------------------------------
! --- hybrid grid generator, single j-row (part A).
! --- --------------------------------------------

  character(len=256) :: mesg  ! A string for output messages
  integer :: fixlay         ! deepest fixed coordinate layer
  real :: qhrlx( GV%ke+1)   ! relaxation coefficient per timestep [nondim]
  real :: dp0ij( GV%ke)     ! minimum layer thickness [H ~> m or kg m-2]
  real :: dp0cum(GV%ke+1)   ! minimum interface depth [H ~> m or kg m-2]

  real :: Rcv_tgt(GV%ke)    ! Target potential density [R ~> kg m-3]
  real :: temp(GV%ke)       ! A column of potential temperature [degC]
  real :: saln(GV%ke)       ! A column of salinity [ppt]
  real :: Rcv(GV%ke)        ! A column of coordinate potential density [R ~> kg m-3]
  real :: h_col(GV%ke)      ! A column of layer thicknesses [H ~> m or kg m-2]
  real :: p_col(GV%ke)      ! A column of reference pressures [R L2 T-2 ~> Pa]
  real :: tracer(GV%ke,max(ntr,1)) ! Columns of each tracer [Conc]
  real :: h_tot             ! Total thickness of the water column [H ~> m or kg m-2]
  real :: nominalDepth      ! Depth of ocean bottom (positive downward) [H ~> m or kg m-2]
  real :: h_thin            ! A negligibly small thickness to identify essentially
                            ! vanished layers [H ~> m or kg m-2]
  real :: dilate            ! A factor by which to dilate the target positions from z to z* [nondim]

  real :: Th_tot_in, Th_tot_out ! Column integrated temperature [degC H ~> degC m or degC kg m-2]
  real :: Sh_tot_in, Sh_tot_out ! Column integrated salinity [ppt H ~> ppt m or ppt kg m-2]
  real :: Trh_tot_in(max(ntr,1))  ! Initial column integrated tracer amounts [conc H ~> conc m or conc kg m-2]
  real :: Trh_tot_out(max(ntr,1)) ! Final column integrated tracer amounts [conc H ~> conc m or conc kg m-2]

  logical :: debug_conservation ! If true, test for non-conservation.
  logical :: terrain_following  ! True if this column is terrain following.
  integer :: trcflg(max(ntr,1)) ! Hycom tracer type flag for each tracer
  integer :: i, j, k, nk, m

  nk = GV%ke

  ! Set all tracers to be passive.  Setting this to 2 treats a tracer like temperature.
  trcflg(:) = 3

  h_thin = 1e-6*GV%m_to_H
  debug_conservation = .false. !  Set this to true for debugging

  p_col(:) = CS%ref_pressure

  do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1 ; if (G%mask2dT(i,j)>0.) then

    h_tot = 0.0
    do k=1,nk
      ! Rcv_tgt(k) = theta(i,j,k)  ! If a 3-d target density were set up in theta, use that here.
      Rcv_tgt(k) = CS%target_density(k)  ! MOM6 does not yet support 3-d target densities.
      h_col(k) = h(i,j,k)
      h_tot = h_tot + h_col(k)
      temp(k) = tv%T(i,j,k)
      saln(k) = tv%S(i,j,k)
    enddo

    ! This sets the potential density from T and S.
    call calculate_density(temp, saln, p_col, Rcv, tv%eqn_of_state)

    do m=1,ntr ; do k=1,nk
      tracer(k,m) = Reg%Tr(m)%t(i,j,k)
    enddo ; enddo

    ! Store original amounts to test for conservation of temperature, salinity, and tracers.
    if (debug_conservation) then
      Th_tot_in = 0.0 ; Sh_tot_in = 0.0 ; Trh_tot_in(:) = 0.0
      do k=1,nk
        Sh_tot_in = Sh_tot_in + h_col(k)*saln(k)
        Th_tot_in = Th_tot_in + h_col(k)*temp(k)
      enddo
      do m=1,ntr ; do k=1,nk
        Trh_tot_in(m) = Trh_tot_in(m) + h_col(k)*tracer(k,m)
      enddo ; enddo
    endif

    ! The following block of code is used to trigger z* stretching of the targets heights.
    nominalDepth = (G%bathyT(i,j)+G%Z_ref)*GV%Z_to_H
    if (h_tot <= CS%min_dilate*nominalDepth) then
      dilate = CS%min_dilate
    elseif (h_tot >= CS%max_dilate*nominalDepth) then
      dilate = CS%max_dilate
    else
      dilate = h_tot / nominalDepth
    endif

    terrain_following = (h_tot < dilate*CS%dpns) .and. (CS%dpns >= CS%dsns)

    ! Convert the regridding parameters into specific constraints for this column.
    call hybgen_column_init(nk, CS%nsigma, CS%dp0k, CS%ds0k, CS%dp00i, &
                            CS%topiso_const, CS%qhybrlx, CS%dpns, CS%dsns, h_tot, dilate, &
                            h_col, fixlay, qhrlx, dp0ij, dp0cum)

    ! Do any unmixing of the column that is needed to move the layer properties toward their targets.
    call hybgen_column_unmix(CS, nk, Rcv_tgt, temp, saln, Rcv, tv%eqn_of_state, &
                             ntr, tracer, trcflg, fixlay, qhrlx, h_col, &
                             terrain_following, h_thin)

    ! Store the output from hybgen_unmix in the 3-d arrays.
    do k=1,nk
      h(i,j,k) = h_col(k)
    enddo
    ! Note that temperature and salinity are among the tracers unmixed here.
    do m=1,ntr ; do k=1,nk
      Reg%Tr(m)%t(i,j,k) = tracer(k,m)
    enddo ; enddo
    ! However, temperature and salinity may have been treated differently from other tracers.
    do k=1,nk
      tv%T(i,j,k) = temp(k)
      tv%S(i,j,k) = saln(k)
    enddo

    ! Test for conservation of temperature, salinity, and tracers.
    if (debug_conservation) then
      Th_tot_out = 0.0 ; Sh_tot_out = 0.0 ; Trh_tot_out(:) = 0.0
      do k=1,nk
        Sh_tot_out = Sh_tot_out + h_col(k)*saln(k)
        Th_tot_out = Th_tot_out + h_col(k)*temp(k)
      enddo
      do m=1,ntr ; do k=1,nk
        Trh_tot_out(m) = Trh_tot_out(m) + h_col(k)*tracer(k,m)
      enddo ; enddo
      if (abs(Sh_tot_in - Sh_tot_out) > 1.e-15*(abs(Sh_tot_in) + abs(Sh_tot_out))) then
        write(mesg, '("i,j=",2i8,"Sh_tot = ",2es17.8," err = ",es13.4)') &
              i, j, Sh_tot_in, Sh_tot_out, (Sh_tot_in - Sh_tot_out)
        call MOM_error(FATAL, "Mismatched column salinity in hybgen_unmix: "//trim(mesg))
      endif
      if (abs(Th_tot_in - Th_tot_out) > 1.e-10*(abs(Th_tot_in) + abs(Th_tot_out))) then
        write(mesg, '("i,j=",2i8,"Th_tot = ",2es17.8," err = ",es13.4)') &
              i, j, Th_tot_in, Th_tot_out, (Th_tot_in - Th_tot_out)
        call MOM_error(FATAL, "Mismatched column temperature in hybgen_unmix: "//trim(mesg))
      endif
      do m=1,ntr
        if (abs(Trh_tot_in(m) - Trh_tot_out(m)) > 1.e-10*(abs(Trh_tot_in(m)) + abs(Trh_tot_out(m)))) then
          write(mesg, '("i,j=",2i8,"Trh_tot(",i2,") = ",2es17.8," err = ",es13.4)') &
                i, j, m, Trh_tot_in(m), Trh_tot_out(m), (Trh_tot_in(m) - Trh_tot_out(m))
          call MOM_error(FATAL, "Mismatched column tracer in hybgen_unmix: "//trim(mesg))
        endif
      enddo
    endif
  endif ; enddo ; enddo !i & j.

end subroutine hybgen_unmix


!> Unmix the properties in the lowest layer if it is too light.
subroutine hybgen_column_unmix(CS, nk, Rcv_tgt, temp, saln, Rcv, eqn_of_state, &
                               ntr, tracer, trcflg, fixlay, qhrlx, h_col, &
                               terrain_following, h_thin)
  type(hybgen_unmix_CS), intent(in) :: CS  !< hybgen unmixing control structure
  integer,        intent(in)    :: nk           !< The number of layers
  integer,        intent(in)    :: fixlay       !< deepest fixed coordinate layer
  real,           intent(in)    :: qhrlx(nk+1)  !< Relaxation fraction per timestep [nondim], < 1.
  real,           intent(in)    :: Rcv_tgt(nk)  !< Target potential density [R ~> kg m-3]
  real,           intent(inout) :: temp(nk)     !< A column of potential temperature [degC]
  real,           intent(inout) :: saln(nk)     !< A column of salinity [ppt]
  real,           intent(inout) :: Rcv(nk)      !< Coordinate potential density [R ~> kg m-3]
  type(EOS_type), intent(in)    :: eqn_of_state !< Equation of state structure
  integer,        intent(in)    :: ntr          !< The number of registered passive tracers
  real,           intent(inout) :: tracer(nk, max(ntr,1)) !< Columns of the passive tracers [Conc]
  integer,        intent(in)    :: trcflg(max(ntr,1)) !< Hycom tracer type flag for each tracer
  real,           intent(inout) :: h_col(nk+1)  !< Layer thicknesses [H ~> m or kg m-2]
  logical,        intent(in)    :: terrain_following !< True if this column is terrain following
  real,           intent(in)    :: h_thin       !< A negligibly small thickness to identify
                                                !! essentially vanished layers [H ~> m or kg m-2]

!
! --- ------------------------------------------------------------------
! --- hybrid grid generator, single column - ummix lowest massive layer.
! --- ------------------------------------------------------------------
!
  ! Local variables
  real :: h_hat       ! A portion of a layer to move across an interface [H ~> m or kg m-2]
  real :: delt, deltm ! Temperature differences between successive layers [degC]
  real :: dels, delsm ! Salinity differences between successive layers [ppt]
  real :: abs_dRdT    ! The absolute value of the derivative of the coordinate density
                      ! with temperature [R degC-1 ~> kg m-3 degC-1]
  real :: abs_dRdS    ! The absolute value of the derivative of the coordinate density
                      ! with salinity [R ppt-1 ~> kg m-3 ppt-1]
  real :: q, qts      ! Nondimensional fractions in the range of 0 to 1 [nondim]
  real :: frac_dts    ! The fraction of the temperature or salinity difference between successive
                      ! layers by which the source layer's property changes by the loss of water
                      ! that matches the destination layers properties via unmixing [nondim].
  real :: qtr         ! The fraction of the water that will come from the layer below,
                      ! used for updating the concentration of passive tracers [nondim]
  real :: swap_T      ! A swap variable for temperature [degC]
  real :: swap_S      ! A swap variable for salinity [ppt]
  real :: swap_R      ! A swap variable for the coordinate potential density [R ~> kg m-3]
  real :: swap_tr     ! A temporary swap variable for the tracers [conc]
  logical, parameter :: lunmix=.true.     ! unmix a too light deepest layer
  integer :: k, ka, kp, kt, m

  ! --- identify the deepest layer kp with significant thickness (> h_thin)
  kp = 2  !minimum allowed value
  do k=nk,3,-1
    if (h_col(k) >= h_thin) then
      kp = k
      exit
    endif
  enddo !k

  k  = kp  !at least 2
  ka = max(k-2,1)  !k might be 2
!
  if ( ((k > fixlay+1) .and. (.not.terrain_following)) .and. & ! layer not fixed depth
       (h_col(k-1) >= h_thin)        .and. & ! layer above not too thin
       (Rcv_tgt(k) > Rcv(k))   .and. & ! layer is lighter than its target
       ((Rcv(k-1) > Rcv(k)) .and. (Rcv(ka) > Rcv(k))) ) then
!
! ---   water in the deepest inflated layer with significant thickness
! ---   (kp) is too light, and it is lighter than the two layers above.
! ---
! ---   this should only occur when relaxing or nudging layer thickness
! ---   and is a bug (bad interaction with tsadvc) even in those cases
! ---
! ---   entrain the entire layer into the one above
!---    note the double negative in T=T-q*(T-T'), equiv. to T=T+q*(T'-T)
    q = h_col(k) / (h_col(k) + h_col(k-1))
    temp(k-1) = temp(k-1) - q*(temp(k-1) - temp(k))
    saln(k-1) = saln(k-1) - q*(saln(k-1) - saln(k))
    call calculate_density(temp(k-1), saln(k-1), CS%ref_pressure, Rcv(k-1), eqn_of_state)

    do m=1,ntr
      tracer(k-1,m) = tracer(k-1,m) - q*(tracer(k-1,m) - tracer(k,m) )
    enddo !m
! ---   entrained the entire layer into the one above, so now kp=kp-1
    h_col(k-1) = h_col(k-1) + h_col(k)
    h_col(k) = 0.0
    kp = k-1
  elseif ( ((k > fixlay+1) .and. (.not.terrain_following)) .and. & ! layer not fixed depth
           (h_col(k-1) >= h_thin) .and. & ! layer above not too thin
           (Rcv_tgt(k) > Rcv(k))  .and. & ! layer is lighter than its target
           (Rcv(k-1) > Rcv(k)) ) then
! ---   water in the deepest inflated layer with significant thickness
! ---   (kp) is too light, and it is lighter than the layer above, but not the layer two above.
! ---
! ---   swap the entire layer with the one above.
    if (h_col(k) <= h_col(k-1)) then
      ! The bottom layer is thinner; swap the entire bottom layer with a portion of the layer above.
      q = h_col(k) / h_col(k-1)  !<=1.0

      swap_T = temp(k-1)
      temp(k-1) = temp(k-1) + q*(temp(k) - temp(k-1))
      temp(k) = swap_T

      swap_S = saln(k-1)
      saln(k-1) = saln(k-1) + q*(saln(k) - saln(k-1))
      saln(k) = swap_S

      Rcv(k) = Rcv(k-1)
      call calculate_density(temp(k-1), saln(k-1), CS%ref_pressure, Rcv(k-1), eqn_of_state)

      do m=1,ntr
        swap_tr = tracer(k-1,m)
        tracer(k-1,m) = tracer(k-1,m) - q * (tracer(k-1,m) - tracer(k,m))
        tracer(k,m) = swap_tr
      enddo !m
    else
      ! The bottom layer is thicker; swap the entire layer above with a portion of the bottom layer.
      q = h_col(k-1) / h_col(k)  !<1.0

      swap_T = temp(k)
      temp(k) = temp(k) + q*(temp(k-1) - temp(k))
      temp(k-1) = swap_T

      swap_S = saln(k)
      saln(k) = saln(k) + q*(saln(k-1) - saln(k))
      saln(k-1) = swap_S

      Rcv(k-1) = Rcv(k)
      call calculate_density(temp(k), saln(k), CS%ref_pressure, Rcv(k), eqn_of_state)

      do m=1,ntr
        swap_tr = tracer(k,m)
        tracer(k,m) = tracer(k,m) + q * (tracer(k-1,m) - tracer(k,m))
        tracer(k-1,m) = swap_tr
      enddo !m
    endif !bottom too light
  endif

  k  = kp  !at least 2
  ka = max(k-2,1)  !k might be 2

  if ( lunmix .and.  & ! usually .true.
       ((k > fixlay+1) .and. (.not.terrain_following)) .and. & ! layer not fixed depth
       (h_col(k-1) >= h_thin)  .and. & ! layer above not too thin
       (Rcv(k) < Rcv_tgt(k))   .and. & ! layer is lighter than its target
       (Rcv(k) > Rcv_tgt(k-1)) .and. & ! layer is denser than the target above
       (abs(Rcv_tgt(k-1) - Rcv(k-1)) < CS%hybiso) .and. & ! layer above is near its target
       (Rcv(k) - Rcv(k-1) > 0.001*(Rcv_tgt(k) - Rcv_tgt(k-1))) ) then
!
! ---   water in the deepest inflated layer with significant thickness (kp) is too
! ---   light but denser than the layer above, with the layer above near-isopycnal
! ---
! ---   split layer into 2 sublayers, one near the desired density
! ---   and one exactly matching the T&S properties of layer k-1.
! ---   To prevent "runaway" T or S, the result satisfies either
! ---     abs(T.k - T.k-1) <= abs(T.k-N - T.k-1) or
! ---     abs(S.k - S.k-1) <= abs(S.k-N - S.k-1) where
! ---     Rcv.k-1 - Rcv.k-N is at least Rcv_tgt(k-1) - Rcv_tgt(k-2)
! ---   It is also limited to a 50% change in layer thickness.

    ka = 1
    do kt=k-2,2,-1
      if ( Rcv(k-1) - Rcv(kt) >= Rcv_tgt(k-1) - Rcv_tgt(k-2) ) then
        ka = kt  !usually k-2
        exit
      endif
    enddo

    delsm = abs(saln(ka) - saln(k-1))
    dels = abs(saln(k-1) - saln(k))
    deltm = abs(temp(ka) - temp(k-1))
    delt = abs(temp(k-1) - temp(k))

    call calculate_density_derivs(temp(k-1), saln(k-1), CS%ref_pressure, abs_dRdT, abs_dRdS, eqn_of_state)
    ! Bound deltm and delsm based on the equation of state and density differences between layers.
    abs_dRdT = abs(abs_dRdT) ; abs_dRdS = abs(abs_dRdS)
    if (abs_dRdT * deltm > Rcv_tgt(k)-Rcv_tgt(k-1)) deltm = (Rcv_tgt(k)-Rcv_tgt(k-1)) / abs_dRdT
    if (abs_dRdS * delsm > Rcv_tgt(k)-Rcv_tgt(k-1)) delsm = (Rcv_tgt(k)-Rcv_tgt(k-1)) / abs_dRdS

    qts = 0.0
    if (qts*dels < min(delsm-dels, dels)) qts = min(delsm-dels, dels) / dels
    if (qts*delt < min(deltm-delt, delt)) qts = min(deltm-delt, delt) / delt

    ! Note that Rcv_tgt(k) > Rcv(k) > Rcv(k-1), and 0 <= qts <= 1.
    ! qhrlx is relaxation coefficient (inverse baroclinic time steps), 0 <= qhrlx <= 1.
    ! This takes the minimum of the two estimates.
    if ((1.0+qts) * (Rcv_tgt(k)-Rcv(k)) < qts * (Rcv_tgt(k)-Rcv(k-1))) then
      q = qhrlx(k) * ((Rcv_tgt(k)-Rcv(k)) / (Rcv_tgt(k)-Rcv(k-1)))
    else
      q = qhrlx(k) * (qts / (1.0+qts)) ! upper sublayer <= 50% of total
    endif
    frac_dts = q / (1.0-q)     ! 0 <= q <= 0.5, so 0 <= frac_dts <= 1

    h_hat = q * h_col(k)
    h_col(k-1) = h_col(k-1) + h_hat
    h_col(k) = h_col(k) - h_hat

    temp(k) = temp(k) + frac_dts * (temp(k) - temp(k-1))
    saln(k) = saln(k) + frac_dts * (saln(k) - saln(k-1))
    call calculate_density(temp(k), saln(k), CS%ref_pressure, Rcv(k), eqn_of_state)

    if ((ntr > 0) .and. (h_hat /= 0.0)) then
      ! qtr is the fraction of the new upper layer from the old lower layer.
      ! The nonconservative original from Hycom: qtr = h_hat / max(h_hat, h_col(k))  !between 0 and 1
      qtr = h_hat / h_col(k-1) ! Between 0 and 1, noting the h_col(k-1) = h_col(k-1) + h_hat above.
      do m=1,ntr
        if (trcflg(m) == 2) then !temperature tracer
          tracer(k,m) = tracer(k,m) + frac_dts * (tracer(k,m) - tracer(k-1,m))
        else !standard tracer - not split into two sub-layers
          tracer(k-1,m) = tracer(k-1,m) + qtr * (tracer(k,m) - tracer(k-1,m))
        endif !trcflg
      enddo !m
    endif !tracers
  endif !too light

!  ! Fill properties of massless or near-massless (thickness < h_thin) layers
!  ! This was in the Hycom verion, but it appears to be unnecessary in MOM6.
!  do k=kp+1,nk
!    ! --- fill thin and massless layers on sea floor with fluid from above
!    Rcv(k) = Rcv(k-1)
!    do m=1,ntr
!      tracer(k,m) = tracer(k-1,m)
!    enddo !m
!    saln(k) = saln(k-1)
!    temp(k) = temp(k-1)
!  enddo !k

end subroutine hybgen_column_unmix

end module MOM_hybgen_unmix
