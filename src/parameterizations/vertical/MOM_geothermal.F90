module MOM_geothermal
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, 2010.                                          *
!*                                                                     *
!*    This file contains the subroutine (geothemal) that implements    *
!*  a geothermal heating at the bottom.  This can be done either in a  *
!*  layered isopycnal mode, in which the heating raises the density of *
!*  the layer to the target density of the layer above, and then moves *
!*  the water into that layer, or in a simple Eulerian mode, in which  *
!*  the bottommost GEOTHERMAL_THICKNESS are heated.  Geothermal heating*
!*  will also provide a buoyant source of bottom TKE that can be used  *
!*  to further mix the near-bottom water.  In cold fresh water lakes   *
!*  where heating increases density, water should be moved into deeper *
!*  layers, but this is not implemented yet.                           *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, buoy, Rml, eaml, ebml, etc.           *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : register_static_field, time_type, diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_io, only : read_data, slasher
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>

public geothermal, geothermal_init, geothermal_end

type, public :: geothermal_CS ; private
  real    :: dRcv_dT_inplace !   The value of dRcv_dT above which (dRcv_dT is
                             ! negative) the water is heated in place instead
                             ! of moving upward between layers, in kg m-3 K-1.
  real, pointer :: geo_heat(:,:) => NULL()  ! The geothermal heat flux, in
                             ! W m-2.
  real    :: geothermal_thick !  The thickness over which geothermal heating is
                             ! applied, in m.
  logical :: apply_geothermal !  If true, geothermal heating will be applied;
                             ! otherwise GEOTHERMAL_SCALE has been set to 0 and
                             ! there is no heat to apply.

  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
end type geothermal_CS

contains

subroutine geothermal(h, tv, dt, ea, eb, G, GV, CS)
  type(ocean_grid_type),                    intent(inout) :: G
  type(verticalGrid_type),                  intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: h
  type(thermo_var_ptrs),                    intent(inout) :: tv
  real,                                     intent(in)    :: dt
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: ea, eb
  type(geothermal_CS),                      pointer       :: CS

!   This subroutine applies geothermal heating, including the movement of water
! between isopycnal layers to match the target densities.  The heating is
! applied to the bottommost layers that occur within ### of the bottom. If
! the partial derivative of the coordinate density with temperature is positive
! or very small, the layers are simply heated in place.  Any heat that can not
! be applied to the ocean is returned (WHERE)?

! Arguments: h - Layer thickness, in m or kg m-2. (Intent in/out)
!                The units of h are referred to as H below.
!  (in/out)  tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      dt - Time increment, in s.
!  (in/out)  ea - The amount of fluid moved downward into a layer; this should
!                 be increased due to mixed layer detrainment, in the same units
!                 as h - usually m or kg m-2 (i.e., H).
!  (in/out)  eb - The amount of fluid moved upward into a layer; this should
!                 be increased due to mixed layer entrainment, in the same units
!                 as h - usually m or kg m-2 (i.e., H).
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 geothermal_init.

!  real :: resid(SZI_(G),SZJ_(G)) !z1l: never been used.

  real, dimension(SZI_(G)) :: &
    heat_rem,  & ! remaining heat (H * degC)
    h_geo_rem, & ! remaining thickness to apply geothermal heating (units of H)
    Rcv_BL,    & ! coordinate density in the deepest variable density layer (kg/m3)
    p_ref        ! coordiante densities reference pressure (Pa)

  real, dimension(2) :: &
    T2, S2, &   ! temp and saln in the present and target layers (degC and ppt)
    dRcv_dT_, & ! partial derivative of coordinate density wrt temp (kg m-3 K-1)
    dRcv_dS_    ! partial derivative of coordinate density wrt saln (kg m-3 ppt-1)

  real :: Angstrom, H_neglect  ! small thicknesses in H
  real :: Rcv           ! coordinate density of present layer (kg m-3)
  real :: Rcv_tgt       ! coordinate density of target layer (kg m-3)
  real :: dRcv          ! difference between Rcv and Rcv_tgt (kg m-3)
  real :: dRcv_dT       ! partial derivative of coordinate density wrt temp
                        ! in the present layer (kg m-3 K-1); usually negative
  real :: h_heated      ! thickness that is being heated (units of H)
  real :: heat_avail    ! heating available for the present layer (units of Kelvin * H)
  real :: heat_in_place ! heating to warm present layer w/o movement between layers (K * H)
  real :: heat_trans    ! heating available to move water from present layer to target layer (K * H)
  real :: heating       ! heating used to move water from present layer to target layer (K * H)
                        ! 0 <= heating <= heat_trans
  real :: h_transfer    ! thickness moved between layers (units of H)
  real :: wt_in_place   ! relative weighting that goes from 0 to 1 (non-dim)
  real :: I_h           ! inverse thickness (units of 1/H)
  real :: dTemp         ! temperature increase in a layer (Kelvin)
  real :: Irho_cp       ! inverse of heat capacity per unit layer volume (units K H m2 J-1)

  logical :: do_i(SZI_(G))
  integer :: i, j, k, is, ie, js, je, nz, k2, i2
  integer :: isj, iej, num_start, num_left, nkmb, k_tgt


  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_geothermal: "//&
         "Module must be initialized before it is used.")
  if (.not.CS%apply_geothermal) return

  nkmb      = GV%nk_rho_varies
  Irho_cp   = 1.0 / (GV%H_to_kg_m2 * tv%C_p)
  Angstrom  = GV%Angstrom
  H_neglect = GV%H_subroundoff
  p_ref(:)  = tv%P_Ref

  if (.not.ASSOCIATED(tv%T)) call MOM_error(FATAL, "MOM geothermal: "//&
      "Geothermal heating can only be applied if T & S are state variables.")

!  do i=is,ie ; do j=js,je
!    resid(i,j) = tv%internal_heat(i,j)
!  enddo ; enddo

!$OMP parallel do default(none) shared(is,ie,js,je,G,GV,CS,dt,Irho_cp,nkmb,tv,    &
!$OMP                                  p_Ref,h,Angstrom,nz,H_neglect,eb)          &
!$OMP                          private(num_start,heat_rem,do_i,h_geo_rem,num_left,&
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

    num_start = 0
    do i=is,ie
      heat_rem(i) = G%mask2dT(i,j) * (CS%geo_heat(i,j) * (dt*Irho_cp))
      do_i(i) = .true. ; if (heat_rem(i) <= 0.0) do_i(i) = .false.
      if (do_i(i)) num_start = num_start + 1
      h_geo_rem(i) = CS%Geothermal_thick * GV%m_to_H
    enddo
    if (num_start == 0) cycle
    num_left = num_start

    ! Find the first and last columns that need to be worked on.
    isj = ie+1 ; do i=is,ie ; if (do_i(i)) then ; isj = i ; exit ; endif ; enddo
    iej = is-1 ; do i=ie,is,-1 ; if (do_i(i)) then ; iej = i ; exit ; endif ; enddo

    if (nkmb > 0) then
      call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_Ref(:), &
                             Rcv_BL(:), isj, iej-isj+1, tv%eqn_of_state)
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
            call calculate_density_derivs(T2(:), S2(:), p_Ref(:), &
                         dRcv_dT_, dRcv_dS_, 1, 2, tv%eqn_of_state)
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

            I_h = 1.0 / (h(i,j,k_tgt) + h_transfer + H_neglect)
            tv%T(i,j,k_tgt) = (h(i,j,k_tgt) * tv%T(i,j,k_tgt) + &
                               (h_transfer * tv%T(i,j,k) + heating)) * I_h
            tv%S(i,j,k_tgt) = (h(i,j,k_tgt) * tv%S(i,j,k_tgt) + &
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
            ! if ((i==isj) .and. (num_left > 0)) then
            !   do i2=isj+1,iej ; if (do_i(i2)) then
            !     isj = i2 ; exit ! Set the new starting value.
            !   endif ; enddo
            ! endif
            ! if ((i==iej) .and. (num_left > 0)) then
            !   do i2=iej-1,isj,-1 ; if (do_i(i2)) then
            !     iej = i2 ; exit ! Set the new ending value.
            !   endif ; enddo
            ! endif
          endif
        endif
      endif ; enddo
      if (num_left <= 0) exit
    enddo ! k-loop

    if (associated(tv%internal_heat)) then ; do i=is,ie
      tv%internal_heat(i,j) = tv%internal_heat(i,j) + GV%H_to_kg_m2 * &
           (G%mask2dT(i,j) * (CS%geo_heat(i,j) * (dt*Irho_cp)) - heat_rem(i))
    enddo ; endif
  enddo ! j-loop

!  do i=is,ie ; do j=js,je
!    resid(i,j) = tv%internal_heat(i,j) - resid(i,j) - GV%H_to_kg_m2 * &
!           (G%mask2dT(i,j) * (CS%geo_heat(i,j) * (dt*Irho_cp)))
!  enddo ; enddo

end subroutine geothermal

subroutine geothermal_init(Time, G, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ctrl), target, intent(inout) :: diag
  type(geothermal_CS),     pointer       :: CS

! Arguments:
!  (in)      Time       - current model time
!  (in)      G          - ocean grid structure
!  (in)      param_file - structure indicating the open file to parse for
!                         model parameter values
!  (in)      diag       - structure used to regulate diagnostic output
!  (in/out)  CS         - pointer pointing to the module control structure

! This include declares and sets the variable "version".
#include "version_variable.h"

  character(len=40)  :: mod = "MOM_geothermal"  ! module name
  character(len=200) :: inputdir, geo_file, filename, geotherm_var
  real :: scale
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
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "GEOTHERMAL_SCALE", scale, &
                 "The constant geothermal heat flux, a rescaling \n"//&
                 "factor for the heat flux read from GEOTHERMAL_FILE, or \n"//&
                 "0 to disable the geothermal heating.", &
                 units="W m-2 or various", default=0.0)
  CS%apply_geothermal = .not.(scale == 0.0)
  if (.not.CS%apply_geothermal) return

  call safe_alloc_ptr(CS%geo_heat, isd, ied, jsd, jed) ; CS%geo_heat(:,:) = 0.0

  call get_param(param_file, mod, "GEOTHERMAL_FILE", geo_file, &
                 "The file from which the geothermal heating is to be \n"//&
                 "read, or blank to use a constant heating rate.", default=" ")
  call get_param(param_file, mod, "GEOTHERMAL_THICKNESS", CS%geothermal_thick, &
                 "The thickness over which to apply geothermal heating.", &
                 units="m", default=0.1)
  call get_param(param_file, mod, "GEOTHERMAL_DRHO_DT_INPLACE", CS%dRcv_dT_inplace, &
                 "The value of drho_dT above which geothermal heating \n"//&
                 "simply heats water in place instead of moving it between \n"//&
                 "isopycnal layers.  This must be negative.", &
                 units="kg m-3 K-1",  default=-0.01)
  if (CS%dRcv_dT_inplace >= 0.0) call MOM_error(FATAL, "geothermal_init: "//&
         "GEOTHERMAL_DRHO_DT_INPLACE must be negative.")

  if (len_trim(geo_file) >= 1) then
    call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    filename = trim(inputdir)//trim(geo_file)
    call log_param(param_file, mod, "INPUTDIR/GEOTHERMAL_FILE", filename)
    call get_param(param_file, mod, "GEOTHERMAL_VARNAME", geotherm_var, &
                 "The name of the geothermal heating variable in \n"//&
                 "GEOTHERMAL_FILE.", default="geo_heat")
    call read_data(filename, trim(geotherm_var), CS%geo_heat, &
                   domain=G%Domain%mpp_domain)
    do j=jsd,jed ; do i=isd,ied
      CS%geo_heat(i,j) = (G%mask2dT(i,j) * scale) * CS%geo_heat(i,j)
    enddo ; enddo
  else
    do j=jsd,jed ; do i=isd,ied
      CS%geo_heat(i,j) = G%mask2dT(i,j) * scale
    enddo ; enddo
  endif

  ! post the static geothermal heating field
  id = register_static_field('ocean_model', 'geo_heat', diag%axesT1,   &
        'Geothermal heat flux into ocean', 'W m-2',                    &
        cmor_field_name='hfgeou', cmor_units='W m-2',                  &
        cmor_standard_name='upward_geothermal_heat_flux_at_sea_floor', &
        cmor_long_name='Upward geothermal heat flux at sea floor')
  if (id > 0) call post_data(id, CS%geo_heat, diag, .true.)

end subroutine geothermal_init

subroutine geothermal_end(CS)
  type(geothermal_CS), pointer :: CS

  if (associated(CS%geo_heat)) deallocate(CS%geo_heat)
  if (associated(CS)) deallocate(CS)
end subroutine geothermal_end

end module MOM_geothermal
