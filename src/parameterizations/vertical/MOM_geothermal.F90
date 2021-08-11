!> Implemented geothermal heating at the ocean bottom.
module MOM_geothermal

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : register_static_field, time_type, diag_ctrl
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_io,            only : MOM_read_data, slasher
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type, get_thickness_units
use MOM_EOS,           only : calculate_density, calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>

public geothermal_entraining, geothermal_in_place, geothermal_init, geothermal_end

!> Control structure for geothermal heating
type, public :: geothermal_CS ; private
  real    :: dRcv_dT_inplace  !< The value of dRcv_dT above which (dRcv_dT is negative) the
                              !! water is heated in place instead of moving upward between
                              !! layers in non-ALE layered mode [R degC-1 ~> kg m-3 degC-1]
  real, allocatable, dimension(:,:) :: geo_heat !< The geothermal heat flux [J m-2 T-1 ~> W m-2]
  real    :: geothermal_thick !< The thickness over which geothermal heating is
                              !! applied [H ~> m or kg m-2]
  logical :: apply_geothermal !< If true, geothermal heating will be applied.  This is false if
                              !! GEOTHERMAL_SCALE is 0 and there is no heat to apply.

  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the timing
                                             !! timing of diagnostic output
  integer :: id_internal_heat_heat_tendency = -1  !< ID for diagnostic of heat tendency
  integer :: id_internal_heat_temp_tendency = -1  !< ID for diagnostic of temperature tendency
  integer :: id_internal_heat_h_tendency = -1     !< ID for diagnostic of thickness tendency

end type geothermal_CS

contains

!> Applies geothermal heating, including the movement of water
!! between isopycnal layers to match the target densities.  The heating is
!! applied to the bottommost layers that occur within GEOTHERMAL_THICKNESS of the bottom. If
!! the partial derivative of the coordinate density with temperature is positive
!! or very small, the layers are simply heated in place.  Any heat that can not
!! be applied to the ocean is returned (WHERE)?
subroutine geothermal_entraining(h, tv, dt, ea, eb, G, GV, US, CS, halo)
  type(ocean_grid_type),                     intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(inout) :: tv !< A structure containing pointers
                                                                 !! to any available thermodynamic fields.
  real,                                      intent(in)    :: dt !< Time increment [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: ea !< The amount of fluid moved
                                                                 !! downward into a layer; this
                                                                 !! should be increased due to mixed
                                                                 !! layer detrainment [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: eb !< The amount of fluid moved upward
                                                                 !! into a layer; this should be
                                                                 !! increased due to mixed layer
                                                                 !! entrainment [H ~> m or kg m-2].
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  type(geothermal_CS),                       pointer       :: CS !< The control structure returned by
                                                                 !! a previous call to
                                                                 !! geothermal_init.
  integer,                         optional, intent(in)    :: halo !< Halo width over which to work
  ! Local variables
  real, dimension(SZI_(G)) :: &
    heat_rem,  & ! remaining heat [H degC ~> m degC or kg degC m-2]
    h_geo_rem, & ! remaining thickness to apply geothermal heating [H ~> m or kg m-2]
    Rcv_BL,    & ! coordinate density in the deepest variable density layer [R ~> kg m-3]
    p_ref        ! coordinate densities reference pressure [R L2 T-2 ~> Pa]

  real, dimension(2) :: &
    T2, S2, &   ! temp and saln in the present and target layers [degC] and [ppt]
    dRcv_dT_, & ! partial derivative of coordinate density wrt temp [R degC-1 ~> kg m-3 degC-1]
    dRcv_dS_    ! partial derivative of coordinate density wrt saln [R ppt-1 ~> kg m-3 ppt-1]

  real :: Angstrom, H_neglect  ! small thicknesses [H ~> m or kg m-2]
  real :: Rcv           ! coordinate density of present layer [R ~> kg m-3]
  real :: Rcv_tgt       ! coordinate density of target layer [R ~> kg m-3]
  real :: dRcv          ! difference between Rcv and Rcv_tgt [R ~> kg m-3]
  real :: dRcv_dT       ! partial derivative of coordinate density wrt temp
                        ! in the present layer [R degC-1 ~> kg m-3 degC-1]; usually negative
  real :: h_heated      ! thickness that is being heated [H ~> m or kg m-2]
  real :: heat_avail    ! heating available for the present layer [degC H ~> degC m or degC kg m-2]
  real :: heat_in_place ! heating to warm present layer w/o movement between layers
                        ! [degC H ~> degC m or degC kg m-2]
  real :: heat_trans    ! heating available to move water from present layer to target
                        ! layer [degC H ~> degC m or degC kg m-2]
  real :: heating       ! heating used to move water from present layer to target layer
                        ! [degC H ~> degC m or degC kg m-2]
                        ! 0 <= heating <= heat_trans
  real :: h_transfer    ! thickness moved between layers [H ~> m or kg m-2]
  real :: wt_in_place   ! relative weighting that goes from 0 to 1 [nondim]
  real :: I_h           ! inverse thickness [H-1 ~> m-1 or m2 kg-1]
  real :: dTemp         ! temperature increase in a layer [degC]
  real :: Irho_cp       ! inverse of heat capacity per unit layer volume
                        ! [degC H Q-1 R-1 Z-1 ~> degC m3 J-1 or degC kg J-1]

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    T_old, & ! Temperature of each layer before any heat is added, for diagnostics [degC]
    h_old, & ! Thickness of each layer before any heat is added, for diagnostics [H ~> m or kg m-2]
    work_3d ! Scratch variable used to calculate changes due to geothermal
  real :: Idt           ! inverse of the timestep [T-1 ~> s-1]

  logical :: do_i(SZI_(G))
  logical :: compute_h_old, compute_T_old
  integer :: i, j, k, is, ie, js, je, nz, k2, i2
  integer :: isj, iej, num_left, nkmb, k_tgt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  endif

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_geothermal: "//&
         "Module must be initialized before it is used.")
  if (.not.CS%apply_geothermal) return

  nkmb      = GV%nk_rho_varies
  Irho_cp   = 1.0 / (GV%H_to_RZ * tv%C_p)
  Angstrom  = GV%Angstrom_H
  H_neglect = GV%H_subroundoff
  p_ref(:)  = tv%P_Ref
  Idt       = 1.0 / dt

  if (.not.associated(tv%T)) call MOM_error(FATAL, "MOM geothermal_entraining: "//&
      "Geothermal heating can only be applied if T & S are state variables.")

!  do j=js,je ; do i=is,ie
!    resid(i,j) = tv%internal_heat(i,j)
!  enddo ; enddo

  ! Conditionals for tracking diagnostic depdendencies
  compute_h_old = CS%id_internal_heat_h_tendency > 0 &
                  .or. CS%id_internal_heat_heat_tendency > 0 &
                  .or. CS%id_internal_heat_temp_tendency > 0

  compute_T_old = CS%id_internal_heat_heat_tendency > 0 &
                  .or. CS%id_internal_heat_temp_tendency > 0

  if (CS%id_internal_heat_heat_tendency > 0) work_3d(:,:,:) = 0.0

  if (compute_h_old .or. compute_T_old) then ; do k=1,nz ; do j=js,je ; do i=is,ie
    ! Save temperature and thickness before any changes are made (for diagnostics)
    h_old(i,j,k) = h(i,j,k)
    T_old(i,j,k) = tv%T(i,j,k)
  enddo ; enddo ; enddo ; endif

!$OMP parallel do default(none) shared(is,ie,js,je,G,GV,US,CS,dt,Irho_cp,nkmb,tv, &
!$OMP                                  p_Ref,h,Angstrom,nz,H_neglect,eb,          &
!$OMP                                  h_old,T_old,work_3d,Idt)                   &
!$OMP                          private(heat_rem,do_i,h_geo_rem,num_left,          &
!$OMP                                  isj,iej,Rcv_BL,h_heated,heat_avail,k_tgt,  &
!$OMP                                  Rcv_tgt,Rcv,dRcv_dT,T2,S2,dRcv_dT_,        &
!$OMP                                  dRcv_dS_,heat_in_place,heat_trans,         &
!$OMP                                  wt_in_place,dTemp,dRcv,h_transfer,heating, &
!$OMP                                  I_h)

  do j=js,je
    ! 1. Only work on columns that are being heated.
    ! 2. Find the deepest layer with any mass.
    ! 3. Find the partial derivative of locally referenced potential density
    !  and coordinate density with temperature, and the density of the layer
    !  and the layer above.
    ! 4. Heat a portion of the bottommost layer until it matches the target
    !    density of the layer above, and move it.
    ! 4a. In the case of variable density layers, heat but do not move.
    ! 5. If there is still heat left over, repeat for the next layer up.
    ! This subroutine updates thickness, T & S, and increments eb accordingly.

    ! 6. If there is not enough mass in the ocean, pass some of the heat up
    !    from the ocean via the frazil field?

    num_left = 0
    do i=is,ie
      heat_rem(i) = G%mask2dT(i,j) * (CS%geo_heat(i,j) * (dt*Irho_cp))
      do_i(i) = .true. ; if (heat_rem(i) <= 0.0) do_i(i) = .false.
      if (do_i(i)) num_left = num_left + 1
      h_geo_rem(i) = CS%Geothermal_thick
    enddo
    if (num_left == 0) cycle

    ! Find the first and last columns that need to be worked on.
    isj = ie+1 ; do i=is,ie ; if (do_i(i)) then ; isj = i ; exit ; endif ; enddo
    iej = is-1 ; do i=ie,is,-1 ; if (do_i(i)) then ; iej = i ; exit ; endif ; enddo

    if (nkmb > 0) then
      call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_Ref(:), Rcv_BL(:), &
                             tv%eqn_of_state, (/isj-(G%isd-1),iej-(G%isd-1)/) )
    else
      Rcv_BL(:) = -1.0
    endif

    do k=nz,1,-1
      do i=isj,iej ; if (do_i(i)) then

        if (h(i,j,k) > Angstrom) then
          if ((h(i,j,k)-Angstrom) >= h_geo_rem(i)) then
            h_heated = h_geo_rem(i)
            heat_avail = heat_rem(i)
            h_geo_rem(i) = 0.0
          else
            h_heated = (h(i,j,k)-Angstrom)
            heat_avail = heat_rem(i) * (h_heated / &
                                        (h_geo_rem(i) + H_neglect))
            h_geo_rem(i) = h_geo_rem(i) - h_heated
          endif

          if (k<=nkmb .or. nkmb<=0) then
            ! Simply heat the layer; convective adjustment occurs later
            ! if necessary.
            k_tgt = k
          elseif ((k==nkmb+1) .or. (GV%Rlay(k-1) < Rcv_BL(i))) then
            ! Add enough heat to match the lowest buffer layer density.
            k_tgt = nkmb
            Rcv_tgt = Rcv_BL(i)
          else
            ! Add enough heat to match the target density of layer k-1.
            k_tgt = k-1
            Rcv_tgt = GV%Rlay(k-1)
          endif

          if (k<=nkmb .or. nkmb<=0) then
            Rcv = 0.0 ; dRcv_dT = 0.0 ! Is this OK?
          else
            call calculate_density(tv%T(i,j,k), tv%S(i,j,k), tv%P_Ref, &
                         Rcv, tv%eqn_of_state)
            T2(1) = tv%T(i,j,k) ; S2(1) = tv%S(i,j,k)
            T2(2) = tv%T(i,j,k_tgt) ; S2(2) = tv%S(i,j,k_tgt)
            call calculate_density_derivs(T2(:), S2(:), p_Ref(:), dRcv_dT_, dRcv_dS_, &
                         tv%eqn_of_state, (/1,2/) )
            dRcv_dT = 0.5*(dRcv_dT_(1) + dRcv_dT_(2))
          endif

          if ((dRcv_dT >= 0.0) .or. (k<=nkmb .or. nkmb<=0)) then
            ! This applies to variable density layers.
            heat_in_place = heat_avail
            heat_trans = 0.0
          elseif (dRcv_dT <= CS%dRcv_dT_inplace) then
            ! This is the option that usually applies in isopycnal coordinates.
            heat_in_place = min(heat_avail, max(0.0, h(i,j,k) * &
                                            ((GV%Rlay(k)-Rcv) / dRcv_dT)))
            heat_trans = heat_avail - heat_in_place
          else
            ! wt_in_place should go from 0 to 1.
            wt_in_place = (CS%dRcv_dT_inplace - dRcv_dT) / CS%dRcv_dT_inplace
            heat_in_place = max(wt_in_place*heat_avail, &
                                h(i,j,k) * ((GV%Rlay(k)-Rcv) / dRcv_dT) )
            heat_trans = heat_avail - heat_in_place
          endif

          if (heat_in_place > 0.0) then
            ! This applies to variable density layers. In isopycnal coordinates
            ! this only arises for relatively fresh water near the freezing
            ! point, in which case heating in place will eventually cause things
            ! to sort themselves out, if only because the water will warm to
            ! the temperature of maximum density.
            dTemp = heat_in_place / (h(i,j,k) + H_neglect)
            tv%T(i,j,k) = tv%T(i,j,k) + dTemp
            heat_rem(i) = heat_rem(i) - heat_in_place
            Rcv = Rcv + dRcv_dT * dTemp
          endif

          if (heat_trans > 0.0) then
            ! The second expression might never be used, but will avoid
            ! division by 0.
            dRcv = max(Rcv - Rcv_tgt, 0.0)

            !   dTemp = -dRcv / dRcv_dT
            !   h_transfer = min(heat_rem(i) / dTemp, h(i,j,k)-Angstrom)
            if ((-dRcv_dT * heat_trans) >= dRcv * (h(i,j,k)-Angstrom)) then
              h_transfer = h(i,j,k) - Angstrom
              heating = (h_transfer * dRcv) / (-dRcv_dT)
              ! Since not all the heat has been applied, return the fraction
              ! of the layer thickness that has not yet been fully heated to
              ! the h_geo_rem.
              h_geo_rem(i) = h_geo_rem(i) + h_heated * &
                      ((heat_avail - (heating + heat_in_place)) / heat_avail)
            else
              h_transfer = (-dRcv_dT * heat_trans) / dRcv
              heating = heat_trans
            endif
            heat_rem(i) = heat_rem(i) - heating

            I_h = 1.0 / ((h(i,j,k_tgt) + H_neglect) + h_transfer)
            tv%T(i,j,k_tgt) = ((h(i,j,k_tgt) + H_neglect) * tv%T(i,j,k_tgt) + &
                               (h_transfer * tv%T(i,j,k) + heating)) * I_h
            tv%S(i,j,k_tgt) = ((h(i,j,k_tgt) + H_neglect) * tv%S(i,j,k_tgt) + &
                               h_transfer * tv%S(i,j,k)) * I_h

            h(i,j,k) = h(i,j,k) - h_transfer
            h(i,j,k_tgt) = h(i,j,k_tgt) + h_transfer
            eb(i,j,k_tgt) = eb(i,j,k_tgt) + h_transfer
            if (k_tgt < k-1) then
              do k2 = k_tgt+1,k-1
                eb(i,j,k2) = eb(i,j,k2) + h_transfer
              enddo
            endif
          endif

          if (heat_rem(i) <= 0.0) then
            do_i(i) = .false. ; num_left = num_left-1
            ! For efficiency, uncomment these?
            ! if ((i==isj) .and. (num_left > 0)) then ; do i2=isj+1,iej ; if (do_i(i2)) then
            !   isj = i2 ; exit ! Set the new starting value.
            ! endif ; enddo ; endif
            ! if ((i==iej) .and. (num_left > 0)) then ; do i2=iej-1,isj,-1 ; if (do_i(i2)) then
            !   iej = i2 ; exit ! Set the new ending value.
            ! endif ; enddo ; endif
          endif
        endif

        ! Calculate heat tendency due to addition and transfer of internal heat
        if (CS%id_internal_heat_heat_tendency > 0) then
          work_3d(i,j,k) = ((GV%H_to_RZ*tv%C_p) * Idt) * (h(i,j,k) * tv%T(i,j,k) - h_old(i,j,k) * T_old(i,j,k))
        endif

      endif ; enddo
      if (num_left <= 0) exit
    enddo ! k-loop

    if (associated(tv%internal_heat)) then ; do i=is,ie
      tv%internal_heat(i,j) = tv%internal_heat(i,j) + GV%H_to_RZ * &
           (G%mask2dT(i,j) * (CS%geo_heat(i,j) * (dt*Irho_cp)) - heat_rem(i))
    enddo ; endif
  enddo ! j-loop

  ! Post diagnostic of 3D tendencies (heat, temperature, and thickness) due to internal heat
  if (CS%id_internal_heat_heat_tendency > 0) then
    call post_data(CS%id_internal_heat_heat_tendency, work_3d, CS%diag, alt_h=h_old)
  endif
  if (CS%id_internal_heat_temp_tendency > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = Idt * (tv%T(i,j,k) - T_old(i,j,k))
    enddo ; enddo ; enddo
    call post_data(CS%id_internal_heat_temp_tendency, work_3d, CS%diag, alt_h=h_old)
  endif
  if (CS%id_internal_heat_h_tendency > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = Idt * (h(i,j,k) - h_old(i,j,k))
    enddo ; enddo ; enddo
    call post_data(CS%id_internal_heat_h_tendency, work_3d, CS%diag, alt_h=h_old)
  endif

!  do j=js,je ; do i=is,ie
!    resid(i,j) = tv%internal_heat(i,j) - resid(i,j) - GV%H_to_RZ * &
!           (G%mask2dT(i,j) * (CS%geo_heat(i,j) * (dt*Irho_cp)))
!  enddo ; enddo

end subroutine geothermal_entraining

!> Applies geothermal heating to the bottommost layers that occur within GEOTHERMAL_THICKNESS of
!! the bottom, by simply heating the water in place.  Any heat that can not be applied to the ocean
!! is returned (WHERE)?
subroutine geothermal_in_place(h, tv, dt, G, GV, US, CS, halo)
  type(ocean_grid_type),                     intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(inout) :: tv !< A structure containing pointers
                                                                 !! to any available thermodynamic fields.
  real,                                      intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  type(geothermal_CS),                       pointer       :: CS !< The control structure returned by
                                                                 !! a previous call to geothermal_init.
  integer,                         optional, intent(in)    :: halo !< Halo width over which to work

  ! Local variables
  real, dimension(SZI_(G)) :: &
    heat_rem,  & ! remaining heat [H degC ~> m degC or kg degC m-2]
    h_geo_rem    ! remaining thickness to apply geothermal heating [H ~> m or kg m-2]

  real :: Angstrom, H_neglect  ! small thicknesses [H ~> m or kg m-2]
  real :: heat_here     ! heating applied to the present layer [degC H ~> degC m or degC kg m-2]
  real :: dTemp         ! temperature increase in a layer [degC]
  real :: Irho_cp       ! inverse of heat capacity per unit layer volume
                        ! [degC H Q-1 R-1 Z-1 ~> degC m3 J-1 or degC kg J-1]

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    dTdt_diag           ! Diagnostic of temperature tendency [degC T-1 ~> degC s-1] which might be
                        ! converted into a layer-integrated heat tendency [Q R Z T-1 ~> W m-2]
  real :: Idt           ! inverse of the timestep [T-1 ~> s-1]
  logical :: do_any     ! True if there is more to be done on the current j-row.
  logical :: calc_diags ! True if diagnostic tendencies are needed.
  integer :: i, j, k, is, ie, js, je, nz, i2, isj, iej

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  endif

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_geothermal: "//&
         "Module must be initialized before it is used.")
  if (.not.CS%apply_geothermal) return

  Irho_cp   = 1.0 / (GV%H_to_RZ * tv%C_p)
  Angstrom  = GV%Angstrom_H
  H_neglect = GV%H_subroundoff
  Idt       = 1.0 / dt

  if (.not.associated(tv%T)) call MOM_error(FATAL, "MOM geothermal_in_place: "//&
      "Geothermal heating can only be applied if T & S are state variables.")

!  do i=is,ie ; do j=js,je
!    resid(i,j) = tv%internal_heat(i,j)
!  enddo ; enddo

  ! Conditionals for tracking diagnostic depdendencies
  calc_diags = (CS%id_internal_heat_heat_tendency > 0) .or. (CS%id_internal_heat_temp_tendency > 0)

  if (calc_diags) dTdt_diag(:,:,:) = 0.0

  !$OMP parallel do default(shared) private(heat_rem,do_any,h_geo_rem,isj,iej,heat_here,dTemp)
  do j=js,je
    ! Only work on columns that are being heated, and heat the near-bottom water.

    ! If there is not enough mass in the ocean, pass some of the heat up
    ! from the ocean via the frazil field?

    do_any = .false.
    do i=is,ie
      heat_rem(i) = G%mask2dT(i,j) * (CS%geo_heat(i,j) * (dt*Irho_cp))
      if (heat_rem(i) > 0.0) do_any = .true.
      h_geo_rem(i) = CS%Geothermal_thick
    enddo
    if (.not.do_any) cycle

    ! Find the first and last columns that need to be worked on.
    isj = ie+1 ; do i=is,ie ; if (heat_rem(i) > 0.0) then ; isj = i ; exit ; endif ; enddo
    iej = is-1 ; do i=ie,is,-1 ; if (heat_rem(i) > 0.0) then ; iej = i ; exit ; endif ; enddo

    do k=nz,1,-1
      do_any = .false.
      do i=isj,iej
        if ((heat_rem(i) > 0.0) .and. (h(i,j,k) > Angstrom)) then
          ! Apply some or all of the remaining heat to this layer.
          ! Convective adjustment occurs outside of this module if necessary.
          if ((h(i,j,k)-Angstrom) >= h_geo_rem(i)) then
            heat_here = heat_rem(i)
            h_geo_rem(i) = 0.0
            heat_rem(i) = 0.0
          else
            heat_here = heat_rem(i) * ((h(i,j,k)-Angstrom) / (h_geo_rem(i) + H_neglect))
            h_geo_rem(i) = h_geo_rem(i) - (h(i,j,k)-Angstrom)
            heat_rem(i) = heat_rem(i) - heat_here
          endif

          dTemp = heat_here / (h(i,j,k) + H_neglect)
          tv%T(i,j,k) = tv%T(i,j,k) + dTemp
          if (calc_diags) dTdt_diag(i,j,k) = dTemp * Idt
        endif

        if (heat_rem(i) > 0.0) do_any= .true.
      enddo

      if (.not.do_any) exit
    enddo ! k-loop

    if (associated(tv%internal_heat)) then ; do i=is,ie
      tv%internal_heat(i,j) = tv%internal_heat(i,j) + GV%H_to_RZ * &
           (G%mask2dT(i,j) * (CS%geo_heat(i,j) * (dt*Irho_cp)) - heat_rem(i))
    enddo ; endif
  enddo ! j-loop

  ! Post diagnostics of 3D tendencies of heat and temperature due to geothermal heat
  if (CS%id_internal_heat_temp_tendency > 0) then
    call post_data(CS%id_internal_heat_temp_tendency, dTdt_diag, CS%diag, alt_h=h)
  endif
  if (CS%id_internal_heat_heat_tendency > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      ! Dangerously reuse dTdt_diag for a related variable with different units, going from
      ! units of [degC T-1 ~> degC s-1] to units of [Q R Z T-1 ~> W m-2]
      dTdt_diag(i,j,k) = (GV%H_to_RZ*tv%C_p) * (h(i,j,k) * dTdt_diag(i,j,k))
    enddo ; enddo ; enddo
    call post_data(CS%id_internal_heat_heat_tendency, dTdt_diag, CS%diag, alt_h=h)
  endif

!  do j=js,je ; do i=is,ie
!    resid(i,j) = tv%internal_heat(i,j) - resid(i,j) - GV%H_to_RZ * &
!           (G%mask2dT(i,j) * (CS%geo_heat(i,j) * (dt*Irho_cp)))
!  enddo ; enddo

end subroutine geothermal_in_place

!> Initialize parameters and allocate memory associated with the geothermal heating module.
subroutine geothermal_init(Time, G, GV, US, param_file, diag, CS, useALEalgorithm)
  type(time_type), target, intent(in)    :: Time !< Current model time.
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< Structure used to regulate diagnostic output.
  type(geothermal_CS),     pointer       :: CS   !< Pointer pointing to the module control
                                                 !! structure.
  logical,       optional, intent(in)    :: useALEalgorithm  !< logical for whether to use ALE remapping

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_geothermal"  ! module name
  character(len=48)  :: thickness_units
  ! Local variables
  character(len=200) :: inputdir, geo_file, filename, geotherm_var
  real :: geo_scale  ! A constant heat flux or dimensionally rescaled geothermal flux scaling factor
                     ! [Q R Z T-1 ~> W m-2] or [Q R Z m2 s J-1 T-1 ~> 1]
  integer :: i, j, isd, ied, jsd, jed, id
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (associated(CS)) then
    call MOM_error(WARNING, "geothermal_init called with an associated"// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag
  CS%Time => Time

  ! write parameters to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "GEOTHERMAL_SCALE", geo_scale, &
                 "The constant geothermal heat flux, a rescaling "//&
                 "factor for the heat flux read from GEOTHERMAL_FILE, or "//&
                 "0 to disable the geothermal heating.", &
                 units="W m-2 or various", default=0.0, scale=US%W_m2_to_QRZ_T)
  CS%apply_geothermal = .not.(geo_scale == 0.0)
  if (.not.CS%apply_geothermal) return

  call safe_alloc_alloc(CS%geo_heat, isd, ied, jsd, jed) ; CS%geo_heat(:,:) = 0.0

  call get_param(param_file, mdl, "GEOTHERMAL_FILE", geo_file, &
                 "The file from which the geothermal heating is to be "//&
                 "read, or blank to use a constant heating rate.", default=" ")
  call get_param(param_file, mdl, "GEOTHERMAL_THICKNESS", CS%geothermal_thick, &
                 "The thickness over which to apply geothermal heating.", &
                 units="m", default=0.1, scale=GV%m_to_H)
  call get_param(param_file, mdl, "GEOTHERMAL_DRHO_DT_INPLACE", CS%dRcv_dT_inplace, &
                 "The value of drho_dT above which geothermal heating "//&
                 "simply heats water in place instead of moving it between "//&
                 "isopycnal layers.  This must be negative.", &
                 units="kg m-3 K-1", scale=US%kg_m3_to_R, default=-0.01, &
                 do_not_log=((GV%nk_rho_varies<=0).or.(GV%nk_rho_varies>=GV%ke)) )
  if (CS%dRcv_dT_inplace >= 0.0) call MOM_error(FATAL, "geothermal_init: "//&
         "GEOTHERMAL_DRHO_DT_INPLACE must be negative.")

  if (len_trim(geo_file) >= 1) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    filename = trim(inputdir)//trim(geo_file)
    call log_param(param_file, mdl, "INPUTDIR/GEOTHERMAL_FILE", filename)
    call get_param(param_file, mdl, "GEOTHERMAL_VARNAME", geotherm_var, &
                 "The name of the geothermal heating variable in GEOTHERMAL_FILE.", &
                 default="geo_heat")
    call MOM_read_data(filename, trim(geotherm_var), CS%geo_heat, G%Domain)
    do j=jsd,jed ; do i=isd,ied
      CS%geo_heat(i,j) = (G%mask2dT(i,j) * geo_scale) * CS%geo_heat(i,j)
    enddo ; enddo
  else
    do j=jsd,jed ; do i=isd,ied
      CS%geo_heat(i,j) = G%mask2dT(i,j) * geo_scale
    enddo ; enddo
  endif
  call pass_var(CS%geo_heat, G%domain)

  thickness_units = get_thickness_units(GV)

  ! post the static geothermal heating field
  id = register_static_field('ocean_model', 'geo_heat', diag%axesT1,   &
        'Geothermal heat flux into ocean', 'W m-2', conversion=US%QRZ_T_to_W_m2, &
        cmor_field_name='hfgeou', &
        cmor_standard_name='upward_geothermal_heat_flux_at_sea_floor', &
        cmor_long_name='Upward geothermal heat flux at sea floor', &
        x_cell_method='mean', y_cell_method='mean', area_cell_method='mean')
  if (id > 0) call post_data(id, CS%geo_heat, diag, .true.)

  ! Diagnostic for tendencies due to internal heat (in 3d)
  CS%id_internal_heat_heat_tendency=register_diag_field('ocean_model', &
        'internal_heat_heat_tendency', diag%axesTL, Time,              &
        'Heat tendency (in 3D) due to internal (geothermal) sources',  &
        'W m-2', conversion=US%QRZ_T_to_W_m2, v_extensive=.true.)
  CS%id_internal_heat_temp_tendency=register_diag_field('ocean_model', &
        'internal_heat_temp_tendency', diag%axesTL, Time,              &
        'Temperature tendency (in 3D) due to internal (geothermal) sources', &
        'degC s-1', conversion=US%s_to_T, v_extensive=.true.)
  if (present(useALEalgorithm)) then ; if (.not.useALEalgorithm) then
    ! Do not offer this diagnostic if heating will be in place.
    CS%id_internal_heat_h_tendency=register_diag_field('ocean_model',    &
        'internal_heat_h_tendency', diag%axesTL, Time,                &
        'Thickness tendency (in 3D) due to internal (geothermal) sources', &
        trim(thickness_units)//' s-1', conversion=GV%H_to_MKS*US%s_to_T, v_extensive=.true.)
  endif ; endif

end subroutine geothermal_init

!> Clean up and deallocate memory associated with the geothermal heating module.
subroutine geothermal_end(CS)
  type(geothermal_CS), intent(inout) :: CS !< Geothermal heating control structure that
                                           !! will be deallocated in this subroutine.
  if (allocated(CS%geo_heat)) deallocate(CS%geo_heat)
end subroutine geothermal_end

!> \namespace mom_geothermal
!!
!! Geothermal heating can be added either in a layered isopycnal mode, in which the heating raises the density
!! of the layer to the target density of the layer above, and then moves the water into that layer, or in a
!! simple Eulerian mode, in which the bottommost GEOTHERMAL_THICKNESS are heated.  Geothermal heating will also
!! provide a buoyant source of bottom TKE that can be used to further mix the near-bottom water. In cold fresh
!! water lakes where heating increases density, water should be moved into deeper layers, but this is not
!! implemented yet.

end module MOM_geothermal
