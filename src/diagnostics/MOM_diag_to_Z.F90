!> Maps tracers and velocities into depth space for output as diagnostic quantities.
!!
!! Currently, a piecewise linear subgrid structure is used for tracers, while velocities can
!! use either piecewise constant or piecewise linear structures.
module MOM_diag_to_Z

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domains,       only : pass_var
use MOM_coms,          only : reproducing_sum
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type, diag_axis_init
use MOM_diag_mediator, only : axes_grp, define_axes_group
use MOM_diag_mediator, only : ocean_register_diag
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : slasher, vardesc, query_vardesc, modify_vardesc
use MOM_spatial_means, only : global_layer_mean
use MOM_variables,     only : p3d, p2d
use MOM_verticalGrid,  only : verticalGrid_type

use netcdf

implicit none ; private

#include <MOM_memory.h>

public calculate_Z_diag_fields
public register_Z_tracer
public MOM_diag_to_Z_init
public calculate_Z_transport
public MOM_diag_to_Z_end
public ocean_register_diag_with_z
public find_overlap
public find_limited_slope
public register_Zint_diag
public calc_Zint_diags

!> The control structure for the MOM_diag_to_Z module
type, public :: diag_to_Z_CS ; private
  ! The following arrays are used to store diagnostics calculated in this
  ! module and unavailable outside of it.

  real, pointer, dimension(:,:,:) :: &
    u_z  => NULL(), &   !< zonal velocity remapped to depth space (m/s)
    v_z  => NULL(), &   !< meridional velocity remapped to depth space (m/s)
    uh_z => NULL(), &   !< zonal transport remapped to depth space (m3/s or kg/s)
    vh_z => NULL()      !< meridional transport remapped to depth space (m3/s or kg/s)

  type(p3d) :: tr_z(MAX_FIELDS_)     !< array of tracers, remapped to depth space
  type(p3d) :: tr_model(MAX_FIELDS_) !< pointers to an array of tracers

  real :: missing_vel             = -1.0e34 !< Missing variable fill values for velocities
  real :: missing_trans           = -1.0e34 !< Missing variable fill values for transports
  real :: missing_tr(MAX_FIELDS_) = -1.0e34 !< Missing variable fill values for tracers
  real :: missing_value           = -1.0e34 !< Missing variable fill values for other diagnostics

  integer :: id_u_z  = -1  !< Diagnostic ID for zonal velocity
  integer :: id_v_z  = -1  !< Diagnostic ID for meridional velocity
  integer :: id_uh_Z = -1  !< Diagnostic ID for zonal transports
  integer :: id_vh_Z = -1  !< Diagnostic ID for meridional transports
  integer :: id_tr(MAX_FIELDS_) = -1  !< Diagnostic IDs for tracers
  integer :: id_tr_xyave(MAX_FIELDS_) = -1  !< Diagnostic IDs for spatially averaged tracers

  integer :: num_tr_used = 0 !< Th enumber of tracers in use.
  integer :: nk_zspace = -1 !< The number of levels in the z-space output

  real, pointer :: Z_int(:) => NULL()  !< interface depths of the z-space file (meter)

  !>@{ Axis groups for z-space diagnostic output
  type(axes_grp) :: axesBz,  axesTz,  axesCuz,  axesCvz
  type(axes_grp) :: axesBzi, axesTzi, axesCuzi, axesCvzi
  type(axes_grp) :: axesZ
  !!@}
  integer, dimension(1) :: axesz_out

  type(diag_ctrl), pointer :: diag => NULL() ! A structure that is used to
                                   ! regulate the timing of diagnostic output.

end type diag_to_Z_CS

integer, parameter :: NO_ZSPACE = -1 !< Flag to enable z-space?

contains

!> Return the global horizontal mean in z-space
function global_z_mean(var,G,CS,tracer)
  type(ocean_grid_type), intent(in)  :: G    !< The ocean's grid structure
  type(diag_to_Z_CS),    pointer     :: CS   !< Control structure returned by
                                             !! previous call to diag_to_Z_init.
  real, dimension(SZI_(G), SZJ_(G), CS%nk_zspace), &
                         intent(in)  :: var    !< An array with the variable to average
  integer,               intent(in)  :: tracer !< The tracer index being worked on
  ! Local variables
  real, dimension(SZI_(G), SZJ_(G), CS%nk_zspace)  :: tmpForSumming, weight
  real, dimension(CS%nk_zspace)                    :: global_z_mean, scalarij, weightij
  real, dimension(CS%nk_zspace)                    :: global_temp_scalar, global_weight_scalar
  real :: valid_point, depth_weight
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  nz = CS%nk_zspace

  ! Initialize local arrays
  tmpForSumming(:,:,:) = 0. ; weight(:,:,:) = 0.

  do k=1,nz ; do j=js,je ; do i=is,ie
    valid_point = 1.0
    ! Weight factor for partial bottom cells
    depth_weight = min( max( (-G%Zd_to_m*G%bathyT(i,j)), CS%Z_int(k+1) ) - CS%Z_int(k), 0.)

    ! Flag the point as invalid if it contains missing data, or is below the bathymetry
    if (var(i,j,k) == CS%missing_tr(tracer)) valid_point = 0.
    if (depth_weight == 0.) valid_point = 0.

    weight(i,j,k) = depth_weight * ( (valid_point * (G%areaT(i,j) * G%mask2dT(i,j))) )

    ! If the point is flagged, set the variable itself to zero to avoid NaNs
    if (valid_point == 0.) then
      tmpForSumming(i,j,k) = 0.0
    else
      tmpForSumming(i,j,k) = var(i,j,k) * weight(i,j,k)
    endif
  enddo ; enddo ; enddo

  global_temp_scalar   = reproducing_sum(tmpForSumming,sums=scalarij)
  global_weight_scalar = reproducing_sum(weight,sums=weightij)

  do k=1, nz
    if (scalarij(k) == 0) then
        global_z_mean(k) = 0.0
    else
        global_z_mean(k) = scalarij(k) / weightij(k)
    endif
  enddo

end function global_z_mean

!> This subroutine maps tracers and velocities into depth space for diagnostics.
subroutine calculate_Z_diag_fields(u, v, h, ssh_in, frac_shelf_h, G, GV, CS)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: u    !< The zonal velocity, in m s-1.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(in)    :: v    !< The meridional velocity, in m s-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: h    !< Layer thicknesses, in H (usually m or kg m-2).
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: ssh_in !< Sea surface height in meters.
  real, dimension(:,:),    pointer       :: frac_shelf_h !< The fraction of the cell area covered by
                                                 !! ice shelf, or unassocatiaed if there is no shelf
  type(diag_to_Z_CS),      pointer       :: CS   !< Control structure returned by a previous call
                                                 !! to diag_to_Z_init.
  ! Local variables
  ! Note the deliberately reversed axes in h_f, u_f, v_f, and tr_f.
  real :: ssh(SZI_(G),SZJ_(G))   ! copy of ssh_in (meter or kg/m2)
  real :: e(SZK_(G)+2)           ! z-star interface heights (meter or kg/m2)
  real :: h_f(SZK_(G)+1,SZI_(G)) ! thicknesses of massive layers (meter or kg/m2)
  real :: u_f(SZK_(G)+1,SZIB_(G))! zonal velocity component in any massive layer
  real :: v_f(SZK_(G)+1,SZI_(G)) ! meridional velocity component in any massive layer

  real :: tr_f(SZK_(G),max(CS%num_tr_used,1),SZI_(G)) ! tracer concentration in massive layers
  integer :: nk_valid(SZIB_(G))  ! number of massive layers in a column

  real :: D_pt(SZIB_(G))        ! bottom depth (meter or kg/m2)
  real :: shelf_depth(SZIB_(G)) ! ice shelf depth (meter or kg/m2)
  real :: htot           ! summed layer thicknesses (meter or kg/m2)
  real :: dilate         ! proportion by which to dilate every layer
  real :: wt(SZK_(G)+1)  ! fractional weight for each layer in the
                         ! range between k_top and k_bot (nondim)
  real :: z1(SZK_(G)+1)  ! z1 and z2 are the depths of the top and bottom
  real :: z2(SZK_(G)+1)  ! limits of the part of a layer that contributes
                         ! to a depth level, relative to the cell center
                         ! and normalized by the cell thickness (nondim)
                         ! Note that -1/2 <= z1 < z2 <= 1/2.
  real :: sl_tr(max(CS%num_tr_used,1)) ! normalized slope of the tracer
                                       ! within the cell, in tracer units
  real :: Angstrom ! A minimal layer thickness, in H.
  real :: slope ! normalized slope of a variable within the cell

  real :: layer_ave(CS%nk_zspace)

  logical :: linear_velocity_profiles, ice_shelf

  integer :: k_top, k_bot, k_bot_prev
  integer :: i, j, k, k2, kz, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nk, m, nkml
  integer :: IsgB, IegB, JsgB, JegB
  is   = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nk = G%ke
  Isq  = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  IsgB = G%IsgB ; IegB = G%IegB ; JsgB = G%JsgB ; JegB = G%JegB
  nkml = max(GV%nkml, 1)
  Angstrom = GV%Angstrom_H
  ssh(:,:) = ssh_in
  linear_velocity_profiles = .true.
  ! Update the halos
  call pass_var(ssh, G%Domain)

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "diagnostic_fields_zstar: Module must be initialized before it is used.")

  ice_shelf = associated(frac_shelf_h)

  ! If no fields are needed, return
  if ((CS%id_u_z <= 0) .and. (CS%id_v_z <= 0) .and. (CS%num_tr_used < 1)) return

  ! zonal velocity component
  if (CS%id_u_z > 0) then

    do kz=1,CS%nk_zspace ; do j=js,je ; do I=Isq,Ieq
      CS%u_z(I,j,kz) = CS%missing_vel
    enddo ; enddo ; enddo


    do j=js,je
      shelf_depth(:) = 0. ! initially all is open ocean
      ! Remove all massless layers.
      do I=Isq,Ieq
        nk_valid(I) = 0
        D_pt(I) = 0.5*G%Zd_to_m*(G%bathyT(i+1,j)+G%bathyT(i,j))
        if (ice_shelf) then
          if (frac_shelf_h(i,j)+frac_shelf_h(i+1,j) > 0.) then ! under shelf
            shelf_depth(I) = abs(0.5*(ssh(i+1,j)+ssh(i,j)))
          endif
        endif
      enddo
      do k=1,nk ; do I=Isq,Ieq
        if ((G%mask2dCu(I,j) > 0.5) .and. (h(i,j,k)+h(i+1,j,k) > 4.0*Angstrom)) then
          nk_valid(I) = nk_valid(I) + 1 ; k2 = nk_valid(I)
          h_f(k2,I) = 0.5*(h(i,j,k)+h(i+1,j,k)) ; u_f(k2,I) = u(I,j,k)
        endif
      enddo ; enddo
      do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.5) then
        ! Add an Angstrom thick layer at the bottom with 0 velocity to impose a
        ! no-slip BBC in the output, if anything but piecewise constant is used.
        nk_valid(I) = nk_valid(I) + 1 ; k2 = nk_valid(I)
        h_f(k2,I) = Angstrom ; u_f(k2,I) = 0.0
        ! GM: D_pt is always slightly larger (by 1E-6 or so) than shelf_depth, so
        ! I consider that the ice shelf is grounded when
        ! shelf_depth(I) + 1.0E-3 > D_pt(i)
        if (ice_shelf .and. shelf_depth(I) + 1.0E-3 > D_pt(i)) nk_valid(I)=0
      endif ; enddo


      do I=Isq,Ieq ; if (nk_valid(I) > 0) then
      ! Calculate the z* interface heights for tracers.
        htot = 0.0 ; do k=1,nk_valid(i) ; htot = htot + h_f(k,i) ; enddo
        dilate = 0.0
        if (htot*GV%H_to_m > 2.0*Angstrom) then
           dilate = MAX((D_pt(i) - shelf_depth(i)),Angstrom)/htot
        endif
        e(nk_valid(i)+1) = -D_pt(i)
        do k=nk_valid(i),1,-1 ; e(K) = e(K+1) + h_f(k,i)*dilate ; enddo

      ! Interpolate each variable into depth space.
        k_bot = 1 ; k_bot_prev = -1
        do kz=1,CS%nk_zspace
          call find_overlap(e, CS%Z_int(kz), CS%Z_int(kz+1), nk_valid(I), &
                            k_bot, k_top, k_bot, wt, z1, z2)
          if (k_top>nk_valid(I)) exit

          !GM if top range that is being map is below the shelf, interpolate
          ! otherwise keep missing_vel
          if (CS%Z_int(kz)<=-shelf_depth(I)) then

            if (linear_velocity_profiles) then
              k = k_top
              if (k /= k_bot_prev) then
                ! Calculate the intra-cell profile.
                slope = 0.0 ! ; curv = 0.0
                if ((k < nk_valid(I)) .and. (k > nkml)) call &
                  find_limited_slope(u_f(:,I), e, slope, k)
              endif
              ! This is the piecewise linear form.
              CS%u_z(I,j,kz) = wt(k) * (u_f(k,I) + 0.5*slope*(z2(k) + z1(k)))
              ! For the piecewise parabolic form add the following...
              !     + C1_3*curv*(z2(k)**2 + z2(k)*z1(k) + z1(k)**2))
              do k=k_top+1,k_bot-1
                CS%u_z(I,j,kz) = CS%u_z(I,j,kz) + wt(k)*u_f(k,I)
              enddo
              if (k_bot > k_top) then ; k = k_bot
                ! Calculate the intra-cell profile.
                slope = 0.0 ! ; curv = 0.0
                if ((k < nk_valid(I)) .and. (k > nkml)) call &
                  find_limited_slope(u_f(:,I), e, slope, k)
                 ! This is the piecewise linear form.
                CS%u_z(I,j,kz) = CS%u_z(I,j,kz) + wt(k) * &
                    (u_f(k,I) + 0.5*slope*(z2(k) + z1(k)))
                ! For the piecewise parabolic form add the following...
                !     + C1_3*curv*(z2(k)**2 + z2(k)*z1(k) + z1(k)**2))
              endif
              k_bot_prev = k_bot
            else ! Use piecewise constant profiles.
              CS%u_z(I,j,kz) = wt(k_top)*u_f(k_top,I)
              do k=k_top+1,k_bot
                CS%u_z(I,j,kz) = CS%u_z(I,j,kz) + wt(k)*u_f(k,I)
              enddo
            endif ! linear profiles
          endif ! below shelf
        enddo ! kz-loop
      endif ; enddo ! I-loop and mask
    enddo ! j-loop

    call post_data(CS%id_u_z, CS%u_z, CS%diag)
  endif

  ! meridional velocity component
  if (CS%id_v_z > 0) then
    do kz=1,CS%nk_zspace ; do J=Jsq,Jeq ; do i=is,ie
      CS%v_z(i,J,kz) = CS%missing_vel
    enddo ; enddo ; enddo

    do J=Jsq,Jeq
      shelf_depth(:) = 0.0 ! initially all is open ocean
      ! Remove all massless layers.
      do i=is,ie
        nk_valid(i) = 0 ; D_pt(i) = 0.5*G%Zd_to_m*(G%bathyT(i,j)+G%bathyT(i,j+1))
        if (ice_shelf) then
          if (frac_shelf_h(i,j)+frac_shelf_h(i,j+1) > 0.) then ! under shelf
            shelf_depth(i) = abs(0.5*(ssh(i,j)+ssh(i,j+1)))
          endif
        endif
      enddo
      do k=1,nk ; do i=is,ie
        if ((G%mask2dCv(i,j) > 0.5) .and. (h(i,j,k)+h(i,j+1,k) > 4.0*Angstrom)) then
          nk_valid(i) = nk_valid(i) + 1 ; k2 = nk_valid(i)
          h_f(k2,i) = 0.5*(h(i,j,k)+h(i,j+1,k)) ; v_f(k2,i) = v(i,j,k)
        endif
      enddo ; enddo
      do i=is,ie ; if (G%mask2dCv(i,j) > 0.5) then
        ! Add an Angstrom thick layer at the bottom with 0 velocity to impose a
        ! no-slip BBC in the output, if anything but piecewise constant is used.
        nk_valid(i) = nk_valid(i) + 1 ; k2 = nk_valid(i)
        h_f(k2,i) = Angstrom ; v_f(k2,i) = 0.0
        if (ice_shelf .and. shelf_depth(i) + 1.0E-3 > D_pt(i)) nk_valid(I)=0
      endif ; enddo

      do i=is,ie ; if (nk_valid(i) > 0) then
      ! Calculate the z* interface heights for tracers.
        htot = 0.0 ; do k=1,nk_valid(i) ; htot = htot + h_f(k,i) ; enddo
        dilate = 0.0
        if (htot > 2.0*Angstrom) then
           dilate = MAX((D_pt(i) - shelf_depth(i)),Angstrom)/htot
        endif
        e(nk_valid(i)+1) = -D_pt(i)
        do k=nk_valid(i),1,-1 ; e(K) = e(K+1) + h_f(k,i)*dilate ; enddo

      ! Interpolate each variable into depth space.
        k_bot = 1 ; k_bot_prev = -1
        do kz=1,CS%nk_zspace
          call find_overlap(e, CS%Z_int(kz), CS%Z_int(kz+1), nk_valid(i), &
                            k_bot, k_top, k_bot, wt, z1, z2)
          if (k_top>nk_valid(i)) exit
          !GM if top range that is being map is below the shelf, interpolate
          ! otherwise keep missing_vel
          if (CS%Z_int(kz)<=-shelf_depth(I)) then
            if (linear_velocity_profiles) then
              k = k_top
              if (k /= k_bot_prev) then
                ! Calculate the intra-cell profile.
                slope = 0.0 ! ; curv = 0.0
                if ((k < nk_valid(i)) .and. (k > nkml)) call &
                  find_limited_slope(v_f(:,i), e, slope, k)
              endif
              ! This is the piecewise linear form.
              CS%v_z(i,J,kz) = wt(k) * (v_f(k,i) + 0.5*slope*(z2(k) + z1(k)))
              ! For the piecewise parabolic form add the following...
              !     + C1_3*curv*(z2(k)**2 + z2(k)*z1(k) + z1(k)**2))
              do k=k_top+1,k_bot-1
                CS%v_z(i,J,kz) = CS%v_z(i,J,kz) + wt(k)*v_f(k,i)
              enddo
              if (k_bot > k_top) then ; k = k_bot
                ! Calculate the intra-cell profile.
                slope = 0.0 ! ; curv = 0.0
                if ((k < nk_valid(i)) .and. (k > nkml)) call &
                  find_limited_slope(v_f(:,i), e, slope, k)
                 ! This is the piecewise linear form.
                CS%v_z(i,J,kz) = CS%v_z(i,J,kz) + wt(k) * &
                    (v_f(k,i) + 0.5*slope*(z2(k) + z1(k)))
                ! For the piecewise parabolic form add the following...
                !     + C1_3*curv*(z2(k)**2 + z2(k)*z1(k) + z1(k)**2))
              endif
              k_bot_prev = k_bot
            else ! Use piecewise constant profiles.
              CS%v_z(i,J,kz) = wt(k_top)*v_f(k_top,i)
              do k=k_top+1,k_bot
                CS%v_z(i,J,kz) = CS%v_z(i,J,kz) + wt(k)*v_f(k,i)
              enddo
            endif ! linear profiles
          endif ! below shelf
          enddo ! kz-loop
      endif ; enddo ! i-loop and mask
    enddo ! J-loop

    call post_data(CS%id_v_z, CS%v_z, CS%diag)
  endif

  ! tracer concentrations
  if (CS%num_tr_used > 0) then

    do m=1,CS%num_tr_used ; do kz=1,CS%nk_zspace ; do j=js,je ; do i=is,ie
      CS%tr_z(m)%p(i,j,kz) = CS%missing_tr(m)
    enddo ; enddo ; enddo ; enddo

    do j=js,je
      shelf_depth(:) = 0.0 ! initially all is open ocean
      ! Remove all massless layers.
      do i=is,ie
        nk_valid(i) = 0 ; D_pt(i) = G%Zd_to_m*G%bathyT(i,j)
        if (ice_shelf) then
          if (frac_shelf_h(i,j) > 0.) then ! under shelf
            shelf_depth(i) = abs(ssh(i,j))
          endif
        endif
      enddo
      do k=1,nk ; do i=is,ie
        if ((G%mask2dT(i,j) > 0.5) .and. (h(i,j,k) > 2.0*Angstrom)) then
          nk_valid(i) = nk_valid(i) + 1 ; k2 = nk_valid(i)
          h_f(k2,i) = h(i,j,k)
          if (ice_shelf .and. shelf_depth(I) + 1.0E-3 > D_pt(i)) nk_valid(I)=0
          do m=1,CS%num_tr_used ; tr_f(k2,m,i) = CS%tr_model(m)%p(i,j,k) ; enddo
        endif
      enddo ; enddo

      do i=is,ie ; if (nk_valid(i) > 0) then
      ! Calculate the z* interface heights for tracers.
        htot = 0.0 ;  do k=1,nk_valid(i) ; htot = htot + h_f(k,i) ; enddo
        dilate = 0.0
        if (htot > 2.0*Angstrom) then
           dilate = MAX((D_pt(i) - shelf_depth(i)),Angstrom)/htot
        endif
        e(nk_valid(i)+1) = -D_pt(i)
        do k=nk_valid(i),1,-1 ; e(K) = e(K+1) + h_f(k,i)*dilate ; enddo

      ! Interpolate each variable into depth space.
        k_bot = 1 ; k_bot_prev = -1
        do kz=1,CS%nk_zspace
          call find_overlap(e, CS%Z_int(kz), CS%Z_int(kz+1), nk_valid(i), &
                            k_bot, k_top, k_bot, wt, z1, z2)
          if (k_top>nk_valid(i)) exit
          if (CS%Z_int(kz)<=-shelf_depth(i)) then
            do m=1,CS%num_tr_used
              k = k_top
              if (k /= k_bot_prev) then
                ! Calculate the intra-cell profile.
                sl_tr(m) = 0.0 ! ; cur_tr(m) = 0.0
                if ((k < nk_valid(i)) .and. (k > nkml)) call &
                  find_limited_slope(tr_f(:,m,i), e, sl_tr(m), k)
              endif
              ! This is the piecewise linear form.
              CS%tr_z(m)%p(i,j,kz) = wt(k) * &
                  (tr_f(k,m,i) + 0.5*sl_tr(m)*(z2(k) + z1(k)))
              ! For the piecewise parabolic form add the following...
              !     + C1_3*cur_tr(m)*(z2(k)**2 + z2(k)*z1(k) + z1(k)**2))
              do k=k_top+1,k_bot-1
                CS%tr_z(m)%p(i,j,kz) = CS%tr_z(m)%p(i,j,kz) + wt(k)*tr_f(k,m,i)
              enddo
              if (k_bot > k_top) then
                k = k_bot
                ! Calculate the intra-cell profile.
                sl_tr(m) = 0.0 ! ; cur_tr(m) = 0.0
                if ((k < nk_valid(i)) .and. (k > nkml)) call &
                  find_limited_slope(tr_f(:,m,i), e, sl_tr(m), k)
                ! This is the piecewise linear form.
                CS%tr_z(m)%p(i,j,kz) = CS%tr_z(m)%p(i,j,kz) + wt(k) * &
                    (tr_f(k,m,i) + 0.5*sl_tr(m)*(z2(k) + z1(k)))
                ! For the piecewise parabolic form add the following...
                !     + C1_3*cur_tr(m)*(z2(k)**2 + z2(k)*z1(k) + z1(k)**2))
              endif
            enddo
            k_bot_prev = k_bot
          endif ! below shelf
        enddo ! kz-loop
      endif ; enddo ! i-loop and mask

    enddo ! j-loop

    do m=1,CS%num_tr_used
      if (CS%id_tr(m) > 0) call post_data(CS%id_tr(m), CS%tr_z(m)%p, CS%diag)
      if (CS%id_tr_xyave(m) > 0) then
        layer_ave = global_z_mean(CS%tr_z(m)%p,G,CS,m)
        call post_data(CS%id_tr_xyave(m), layer_ave, CS%diag)
      endif
    enddo
  endif

end subroutine calculate_Z_diag_fields

!> This subroutine maps horizontal transport into depth space for diagnostic output.
subroutine calculate_Z_transport(uh_int, vh_int, h, dt, G, GV, CS)
  type(ocean_grid_type),                     intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV   !< The ocean's vertical grid
                                                                   !! structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: uh_int !< Time integrated zonal
                                                                   !! transport (m3 or kg).
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: vh_int !< Time integrated meridional
                                                                   !! transport (m3 or kg).
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h    !< Layer thicknesses, in H
                                                                   !! (usually m or kg m-2).
  real,                                      intent(in)    :: dt   !< The time difference in s since
                                                                   !! the last call to this
                                                                   !! subroutine.
  type(diag_to_Z_CS),                        pointer       :: CS   !< Control structure returned by
                                                                   !! previous call to
                                                                   !! diag_to_Z_init.
  ! Local variables
  real, dimension(SZI_(G), SZJ_(G)) :: &
    htot, &        ! total layer thickness (meter or kg/m2)
    dilate         ! nondimensional factor by which to dilate layers to
                   ! convert them into z* space.  (-G%D < z* < 0)

  real, dimension(SZI_(G), max(CS%nk_zspace,1)) :: &
    uh_Z           ! uh_int interpolated into depth space (m3 or kg)
  real, dimension(SZIB_(G), max(CS%nk_zspace,1)) :: &
    vh_Z           ! vh_int interpolated into depth space (m3 or kg)

  real :: h_rem    ! dilated thickness of a layer that has yet to be mapped
                   ! into depth space (meter or kg/m2)
  real :: uh_rem   ! integrated zonal transport of a layer that has yet to be
                   ! mapped into depth space (m3 or kg)
  real :: vh_rem   ! integrated meridional transport of a layer that has yet
                   ! to be mapped into depth space (m3 or kg)
  real :: h_here   ! thickness of a layer that is within the range of the
                   ! current depth level (meter or kg/m2)
  real :: h_above  ! thickness of a layer that is above the current depth
                   ! level (meter or kg.m2)
  real :: uh_here  ! zonal transport of a layer that is attributed to the
                   ! current depth level (m3 or kg)
  real :: vh_here  ! meridional transport of a layer that is attributed to
                   ! the current depth level (m3 or kg)
  real :: Idt      ! inverse of the time step (sec)

  real :: Z_int_above(SZIB_(G)) ! height of the interface atop a layer (meter or kg/m2)

  integer :: kz(SZIB_(G)) ! index of depth level that is being contributed to

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nk, nk_z
  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nk = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "calculate_Z_transport: Module must be initialized before it is used.")
  if ((CS%id_uh_Z <= 0) .and. (CS%id_vh_Z <= 0)) return

  Idt = 1.0 ; if (dt > 0.0) Idt = 1.0 / dt
  nk_z = CS%nk_zspace

  if (nk_z <= 0) return

  ! Determine how much the layers will be dilated in recasting them into z*
  ! coordiantes.  (-G%D < z* < 0).
  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    htot(i,j) = GV%H_subroundoff
  enddo ; enddo
  do k=1,nk ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    htot(i,j) = htot(i,j) + h(i,j,k)
  enddo ; enddo ; enddo
  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    dilate(i,j) = G%Zd_to_m*G%bathyT(i,j) / htot(i,j)
  enddo ; enddo

  ! zonal transport
  if (CS%id_uh_Z > 0) then ; do j=js,je
    do I=Isq,Ieq
      kz(I) = nk_z ; z_int_above(I) = -0.5*G%Zd_to_m*(G%bathyT(i,j)+G%bathyT(i+1,j))
    enddo
    do k=nk_z,1,-1 ; do I=Isq,Ieq
      uh_Z(I,k) = 0.0
      if (CS%Z_int(k) < z_int_above(I)) kz(I) = k-1
    enddo ; enddo
    do k=nk,1,-1 ; do I=Isq,Ieq
      h_rem = 0.5*(dilate(i,j)*h(i,j,k) + dilate(i+1,j)*h(i+1,j,k))
      uh_rem = uh_int(I,j,k)
      z_int_above(I) = z_int_above(I) + h_rem

      do ! Distribute this layer's transport into the depth-levels.
        h_above = z_int_above(I) - CS%Z_int(kz(I))
        if ((kz(I) == 1) .or. (h_above <= 0.0) .or. (h_rem <= 0.0)) then
          ! The entire remaining transport is on this level.
          uh_Z(I,kz(I)) = uh_Z(I,kz(I)) + uh_rem ; exit
        else
          h_here = h_rem - h_above
          uh_here = uh_rem * (h_here / h_rem)

          h_rem = h_rem - h_here ! = h_above
          uh_Z(I,kz(I)) = uh_Z(I,kz(I)) + uh_here
          uh_rem = uh_rem - uh_here
          kz(I) = kz(I) - 1
        endif
      enddo ! End of loop through the target depth-space levels.
    enddo ; enddo
    do k=1,nk_z ; do I=Isq,Ieq
      CS%uh_z(I,j,k) = uh_Z(I,k)*Idt
    enddo ; enddo
  enddo ; endif

  ! meridional transport
  if (CS%id_vh_Z > 0) then ; do J=Jsq,Jeq
    do i=is,ie
      kz(i) = nk_z ; z_int_above(i) = -0.5*G%Zd_to_m*(G%bathyT(i,j)+G%bathyT(i,j+1))
    enddo
    do k=nk_z,1,-1 ; do i=is,ie
      vh_Z(i,k) = 0.0
      if (CS%Z_int(k) < z_int_above(i)) kz(i) = k-1
    enddo ; enddo
    do k=nk,1,-1 ; do i=is,ie
      h_rem = 0.5*(dilate(i,j)*h(i,j,k) + dilate(i,j+1)*h(i,j+1,k))
      vh_rem = vh_int(i,J,k)
      z_int_above(i) = z_int_above(i) + h_rem

      do ! Distribute this layer's transport into the depth-levels.
        h_above = z_int_above(i) - CS%Z_int(kz(i))
        if ((kz(i) == 1) .or. (h_above <= 0.0) .or. (h_rem <= 0.0)) then
          ! The entire remaining transport is on this level.
          vh_Z(i,kz(i)) = vh_Z(i,kz(i)) + vh_rem ; exit
        else
          h_here = h_rem - h_above
          vh_here = vh_rem * (h_here / h_rem)

          h_rem = h_rem - h_here ! = h_above
          vh_Z(i,kz(i)) = vh_Z(i,kz(i)) + vh_here
          vh_rem = vh_rem - vh_here
          kz(i) = kz(i) - 1
        endif
      enddo ! End of loop through the target depth-space levels.
    enddo ; enddo
    do k=1,nk_z ; do i=is,ie
      CS%vh_z(i,J,k) = vh_Z(i,k)*Idt
    enddo ; enddo
  enddo ; endif

  if (CS%id_uh_Z > 0) then
    do k=1,nk_z ; do j=js,je ; do I=Isq,Ieq
      CS%uh_z(i,j,k) = CS%uh_z(i,j,k)*GV%H_to_kg_m2
    enddo ; enddo ; enddo
    call post_data(CS%id_uh_Z, CS%uh_z, CS%diag)
  endif

  if (CS%id_vh_Z > 0) then
    do k=1,nk_z ; do j=Jsq,Jeq ; do I=is,ie
      CS%vh_z(i,j,k) = CS%vh_z(i,j,k)*GV%H_to_kg_m2
    enddo ; enddo ; enddo
    call post_data(CS%id_vh_Z, CS%vh_z, CS%diag)
  endif

end subroutine calculate_Z_transport

!> Determines the layers bounded by interfaces e that overlap
!! with the depth range between Z_top and Z_bot, and the fractional weights
!! of each layer. It also calculates the normalized relative depths of the range
!! of each layer that overlaps that depth range.
subroutine find_overlap(e, Z_top, Z_bot, k_max, k_start, k_top, k_bot, wt, z1, z2)
  real, dimension(:), intent(in)    :: e      !< Column interface heights (meter or kg/m2).
  real,               intent(in)    :: Z_top  !< Top of range being mapped to (meter or kg/m2).
  real,               intent(in)    :: Z_bot  !< Bottom of range being mapped to (meter or kg/m2).
  integer,            intent(in)    :: k_max  !< Number of valid layers.
  integer,            intent(in)    :: k_start !< Layer at which to start searching.
  integer,            intent(inout) :: k_top  !< Indices of top layers that overlap with the depth
                                              !! range.
  integer,            intent(inout) :: k_bot  !< Indices of bottom layers that overlap with the
                                              !! depth range.
  real, dimension(:), intent(out)   :: wt     !< Relative weights of each layer from k_top to k_bot.
  real, dimension(:), intent(out)   :: z1     !< Depth of the top limits of the part of
       !! a layer that contributes to a depth level, relative to the cell center and normalized
       !! by the cell thickness (nondim).  Note that -1/2 <= z1 < z2 <= 1/2.
  real, dimension(:), intent(out)   :: z2     !< Depths of the bottom limit of the part of
       !! a layer that contributes to a depth level, relative to the cell center and normalized
       !! by the cell thickness (nondim).  Note that -1/2 <= z1 < z2 <= 1/2.
  ! Local variables
  real    :: Ih, e_c, tot_wt, I_totwt
  integer :: k

  do k=k_start,k_max ; if (e(K+1)<Z_top) exit ; enddo
  k_top = k
  if (k>k_max) return

  ! Determine the fractional weights of each layer.
  ! Note that by convention, e and Z_int decrease with increasing k.
  if (e(K+1)<=Z_bot) then
    wt(k) = 1.0 ; k_bot = k
    Ih = 0.0 ; if (e(K) /= e(K+1)) Ih = 1.0 / (e(K)-e(K+1))
    e_c = 0.5*(e(K)+e(K+1))
    z1(k) = (e_c - MIN(e(K),Z_top)) * Ih
    z2(k) = (e_c - Z_bot) * Ih
  else
    wt(k) = MIN(e(K),Z_top) - e(K+1) ; tot_wt = wt(k) ! These are always > 0.
    if (e(K) /= e(K+1)) then
      z1(k) = (0.5*(e(K)+e(K+1)) - MIN(e(K), Z_top)) / (e(K)-e(K+1))
    else ; z1(k) = -0.5 ; endif
    z2(k) = 0.5
    k_bot = k_max
    do k=k_top+1,k_max
      if (e(K+1)<=Z_bot) then
        k_bot = k
        wt(k) = e(K) - Z_bot ; z1(k) = -0.5
        if (e(K) /= e(K+1)) then
          z2(k) = (0.5*(e(K)+e(K+1)) - Z_bot) / (e(K)-e(K+1))
        else ; z2(k) = 0.5 ; endif
      else
        wt(k) = e(K) - e(K+1) ; z1(k) = -0.5 ; z2(k) = 0.5
      endif
      tot_wt = tot_wt + wt(k) ! wt(k) is always > 0.
      if (k>=k_bot) exit
    enddo

    I_totwt = 1.0 / tot_wt
    do k=k_top,k_bot ; wt(k) = I_totwt*wt(k) ; enddo
  endif

end subroutine find_overlap

!> This subroutine determines a limited slope for val to be advected with
!! a piecewise limited scheme.
subroutine find_limited_slope(val, e, slope, k)
  real, dimension(:), intent(in)  :: val !< A column of values that are being interpolated.
  real, dimension(:), intent(in)  :: e   !< Column interface heights (meter or kg/m2).
  real,               intent(out) :: slope !< Normalized slope in the intracell distribution of val.
  integer,            intent(in)  :: k   !< Layer whose slope is being determined.
  ! Local variables
  real :: d1, d2

  d1 = 0.5*(e(K-1)-e(K+1)) ; d2 = 0.5*(e(K)-e(K+2))
  if (((val(k)-val(k-1)) * (val(k)-val(k+1)) >= 0.0) .or. (d1*d2 <= 0.0)) then
    slope = 0.0 ! ; curvature = 0.0
  else
    slope = (d1**2*(val(k+1) - val(k)) + d2**2*(val(k) - val(k-1))) * &
            ((e(K) - e(K+1)) / (d1*d2*(d1+d2)))
    ! slope = 0.5*(val(k+1) - val(k-1))
    ! This is S.J. Lin's form of the PLM limiter.
    slope = sign(1.0,slope) * min(abs(slope), &
        2.0*(max(val(k-1),val(k),val(k+1)) - val(k)), &
        2.0*(val(k) - min(val(k-1),val(k),val(k+1))))
    ! curvature = 0.0
  endif

end subroutine find_limited_slope

!> This subroutine calculates interface diagnostics in z-space.
subroutine calc_Zint_diags(h, in_ptrs, ids, num_diags, G, GV, CS)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h    !< Layer thicknesses, in H (usually m or kg m-2).
  type(p3d), dimension(:), intent(in) :: in_ptrs !< Pointers to the diagnostics to be regridded
  integer,   dimension(:), intent(in) :: ids  !< The diagnostic IDs of the diagnostics
  integer,                 intent(in) :: num_diags !< The number of diagnostics to regrid
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure.
  type(diag_to_Z_CS),      pointer    :: CS   !< Control structure returned by
                                              !! previous call to diag_to_Z_init.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),max(CS%nk_zspace+1,1),max(num_diags,1)) :: &
    diag_on_Z  ! diagnostics interpolated to depth space
  real, dimension(SZI_(G),SZK_(G)+1) :: e
  real, dimension(max(num_diags,1),SZI_(G),SZK_(G)+1) :: diag2d

  real, dimension(SZI_(G)) :: &
    htot, &              ! summed layer thicknesses (meter or kg/m2)
    dilate               ! proportion by which to dilate every layer
  real :: wt             ! weighting of the interface above in the
                         ! interpolation to target depths
  integer :: kL(SZI_(G)) ! layer-space index of shallowest interface
                         ! below the target depth

  integer :: i, j, k, k2, kz, is, ie, js, je, nk, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nk = G%ke

  if (num_diags < 1) return
  if (.not.associated(CS)) call MOM_error(FATAL, &
       "calc_Zint_diags: Module must be initialized before it is used.")

  do j=js,je
    ! Calculate the stretched z* interface depths.
    do i=is,ie ; htot(i) = 0.0 ; kL(i) = 1 ; enddo
    do k=1,nk ; do i=is,ie ; htot(i) = htot(i) + h(i,j,k) ; enddo ; enddo
    do i=is,ie
      dilate(i) = 0.0
      if (htot(i)*GV%H_to_m > 0.5) dilate(i) = (G%Zd_to_m*G%bathyT(i,j) - 0.0) / htot(i)
      e(i,nk+1) = -G%Zd_to_m*G%bathyT(i,j)
    enddo
    do k=nk,1,-1 ; do i=is,ie
      e(i,k) = e(i,k+1) + h(i,j,k) * dilate(i)
    enddo ; enddo
    ! e(i,1) should be 0 as a consistency check.

    do k=1,nk+1 ; do i=is,ie ; do m=1,num_diags
      diag2d(m,i,k) = in_ptrs(m)%p(i,j,k)
    enddo ; enddo ; enddo

    do kz=1,CS%nk_zspace+1 ; do i=is,ie
      ! Find the interface below the target Z-file depth, kL.
      if (CS%Z_int(kz) < e(i,nk+1)) then
        kL(i) = nk+2
      else
        do k=kL(i),nk+1 ; if (CS%Z_int(kz) > e(i,k)) exit ; enddo
        kL(i) = k
      endif
      if (kL(i)>1) then
        if (CS%Z_int(kz) > e(i,kL(i)-1)) call MOM_error(FATAL, &
        "calc_Zint_diags: Interface depth mapping is incorrect.")
      endif
      if ((kL(i)>1) .and. (kL(i)<=nk+1)) then
        if (e(i,kL(i)-1) == e(i,kL(i))) call MOM_error(WARNING, &
          "calc_Zint_diags: Interface depths equal.", all_print=.true.)
        if (e(i,kL(i)-1) - e(i,kL(i)) < 0.0) call MOM_error(FATAL, &
          "calc_Zint_diags: Interface depths inverted.")
      endif

      if (kL(i) <= 1) then
        do m=1,num_diags
          diag_on_Z(i,j,kz,m) = diag2d(m,i,1)
        enddo
      elseif (kL(i) > nk+1) then
        do m=1,num_diags
          diag_on_Z(i,j,kz,m) = CS%missing_value
        enddo
      else
        wt = 0.0 ! This probably should not happen?
        if (e(i,kL(i)-1) - e(i,kL(i)) > 0.0) &
          wt = (CS%Z_int(kz) - e(i,kL(i))) / (e(i,kL(i)-1) - e(i,kL(i)))
        if ((wt < 0.0) .or. (wt > 1.0)) call MOM_error(FATAL, &
          "calc_Zint_diags: Bad value of wt found.")
        do m=1,num_diags
          diag_on_Z(i,j,kz,m) = wt * diag2d(m,i,kL(i)-1) + &
                                (1.0-wt) * diag2d(m,i,kL(i))
        enddo
      endif
    enddo ; enddo

  enddo

  do m=1,num_diags
    if (ids(m) > 0) call post_data(ids(m), diag_on_Z(:,:,:,m), CS%diag)
  enddo

end subroutine calc_Zint_diags

!> This subroutine registers a tracer to be output in depth space.
subroutine register_Z_tracer(tr_ptr, name, long_name, units, Time, G, CS, standard_name,   &
     cmor_field_name, cmor_long_name, cmor_units, cmor_standard_name)
  type(ocean_grid_type),      intent(in) :: G      !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                    target,   intent(in) :: tr_ptr !< Tracer for translation to Z-space.
  character(len=*),           intent(in) :: name   !< name for the output tracer.
  character(len=*),           intent(in) :: long_name !< Long name for the output tracer.
  character(len=*),           intent(in) :: units  !< Units of output tracer.
  character(len=*), optional, intent(in) :: standard_name !< The CMOR standard name of this variable.
  type(time_type),            intent(in) :: Time   !< Current model time.
  type(diag_to_Z_CS),         pointer    :: CS     !< Control struct returned by previous
                                                   !! call to diag_to_Z_init.
  character(len=*), optional, intent(in) :: cmor_field_name !< cmor name of a field.
  character(len=*), optional, intent(in) :: cmor_long_name  !< cmor long name of a field.
  character(len=*), optional, intent(in) :: cmor_units      !< cmor units of a field.
  character(len=*), optional, intent(in) :: cmor_standard_name !< cmor standardized name
                                                            !! associated with a field.

  ! Local variables
  character(len=256) :: posted_standard_name
  character(len=256) :: posted_cmor_units
  character(len=256) :: posted_cmor_standard_name
  character(len=256) :: posted_cmor_long_name

  if (CS%nk_zspace<1) return

  if (present(standard_name)) then
    posted_standard_name = standard_name
  else
    posted_standard_name = 'not provided'
  endif

  call register_Z_tracer_low(tr_ptr, name, long_name, units, trim(posted_standard_name), Time, G, CS)

  if (present(cmor_field_name)) then
    ! Fallback values for strings set to "NULL"
    posted_cmor_units         = "not provided"   !
    posted_cmor_standard_name = "not provided"   ! values might be replaced with a CS%missing field?
    posted_cmor_long_name     = "not provided"   !

    ! If attributes are present for MOM variable names, use them first for the register_diag_field
    ! call for CMOR verison of the variable
    posted_cmor_units         = units
    posted_cmor_long_name     = long_name
    posted_cmor_standard_name = posted_standard_name

    ! If specified in the call to register_diag_field, override attributes with the CMOR versions
    if (present(cmor_units)) posted_cmor_units                 = cmor_units
    if (present(cmor_standard_name)) posted_cmor_standard_name = cmor_standard_name
    if (present(cmor_long_name)) posted_cmor_long_name         = cmor_long_name

    call register_Z_tracer_low(tr_ptr, trim(cmor_field_name), trim(posted_cmor_long_name),&
         trim(posted_cmor_units), trim(posted_cmor_standard_name), Time, G, CS)

  endif

end subroutine register_Z_tracer

!> This subroutine registers a tracer to be output in depth space.
subroutine register_Z_tracer_low(tr_ptr, name, long_name, units, standard_name, Time, G, CS)
  type(ocean_grid_type), intent(in) :: G      !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                 target, intent(in) :: tr_ptr !< Tracer for translation to Z-space.
  character(len=*),      intent(in) :: name   !< Name for the output tracer.
  character(len=*),      intent(in) :: long_name !< Long name for output tracer.
  character(len=*),      intent(in) :: units  !< Units of output tracer.
  character(len=*),      intent(in) :: standard_name !< The CMOR standard name of this variable.
  type(time_type),       intent(in) :: Time   !< Current model time.
  type(diag_to_Z_CS),    pointer    :: CS     !< Control struct returned by previous call to
                                              !! diag_to_Z_init.
  ! Local variables
  character(len=256) :: posted_standard_name
  integer :: isd, ied, jsd, jed, nk, m, id_test
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nk = G%ke

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "register_Z_tracer: Module must be initialized before it is used.")

  if (CS%num_tr_used >= MAX_FIELDS_) then
    call MOM_error(WARNING,"MOM_diag_to_Z:  Attempted to register and use "//&
                   "more than MAX_FIELDS_ z-space tracers via register_Z_tracer.")
    return
  endif

  m = CS%num_tr_used + 1

  CS%missing_tr(m) = CS%missing_value ! This could be changed later, if desired.
  if (CS%nk_zspace > 0) then
    CS%id_tr(m) = register_diag_field('ocean_model_zold', name, CS%axesTz, Time,           &
                                      long_name, units, missing_value=CS%missing_tr(m), &
                                      standard_name=standard_name)
    CS%id_tr_xyave(m) = register_diag_field('ocean_model_zold', trim(name)//'_xyave', CS%axesZ, Time, &
                                      long_name, units, missing_value=CS%missing_tr(m), &
                                      standard_name=standard_name)
  else
    id_test = register_diag_field('ocean_model_zold', name, CS%diag%axesT1, Time,      &
                                  long_name, units, missing_value=CS%missing_tr(m), &
                                  standard_name=standard_name)
    if (id_test>0) call MOM_error(WARNING, &
        "MOM_diag_to_Z_init: "//trim(name)// &
        " cannot be output without an appropriate depth-space target file.")
  endif

  if (CS%id_tr(m) <= 0) CS%id_tr(m) = -1
  if (CS%id_tr_xyave(m) <= 0) CS%id_tr_xyave(m) = -1
  if (CS%id_tr(m) > 0 .or. CS%id_tr_xyave(m) > 0) then
    CS%num_tr_used = m
    call safe_alloc_ptr(CS%tr_z(m)%p,isd,ied,jsd,jed,CS%nk_zspace)
    CS%tr_model(m)%p => tr_ptr
  endif

end subroutine register_Z_tracer_low

!> This subroutine sets parameters that control Z-space diagnostic output.
subroutine MOM_diag_to_Z_init(Time, G, GV, param_file, diag, CS)
  type(time_type),         intent(in)    :: Time !< Current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< Struct to regulate diagnostic output.
  type(diag_to_Z_CS),      pointer       :: CS   !< Pointer to point to control structure for
                                                 !! this module, which is allocated and
                                                 !! populated here.
! This include declares and sets the variable "version".
#include "version_variable.h"
  ! Local variables
  character(len=40)  :: mdl = "MOM_diag_to_Z" ! module name
  character(len=200) :: in_dir, zgrid_file    ! strings for directory/file
  character(len=48)  :: flux_units, string
  integer :: z_axis, zint_axis
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nk, id_test
  isd  = G%isd   ; ied = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nk = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "MOM_diag_to_Z_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  if (GV%Boussinesq) then ; flux_units = "m3 s-1"
  else ; flux_units = "kg s-1" ; endif

  CS%diag => diag

  ! Read parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  ! Read in z-space info from a NetCDF file.
  call get_param(param_file, mdl, "Z_OUTPUT_GRID_FILE", zgrid_file, &
                 "The file that specifies the vertical grid for \n"//&
                 "depth-space diagnostics, or blank to disable \n"//&
                 "depth-space output.", default="")

  if (len_trim(zgrid_file) > 0) then
    call get_param(param_file, mdl, "INPUTDIR", in_dir, &
                 "The directory in which input files are found.", default=".")
    in_dir = slasher(in_dir)
    call get_Z_depths(trim(in_dir)//trim(zgrid_file), "zw", CS%Z_int, "zt", &
                      z_axis, zint_axis, CS%nk_zspace)
    call log_param(param_file, mdl, "!INPUTDIR/Z_OUTPUT_GRID_FILE", &
                   trim(in_dir)//trim(zgrid_file))
    call log_param(param_file, mdl, "!NK_ZSPACE (from file)", CS%nk_zspace, &
                 "The number of depth-space levels.  This is determined \n"//&
                 "from the size of the variable zw in the output grid file.", &
                 units="nondim")
  else
    CS%nk_zspace = -1
  endif

  if (CS%nk_zspace > 0) then

    call define_axes_group(diag, (/ diag%axesB1%handles(1),  diag%axesB1%handles(2),  z_axis /),    CS%axesBz)
    call define_axes_group(diag, (/ diag%axesT1%handles(1),  diag%axesT1%handles(2),  z_axis /),    CS%axesTz)
    call define_axes_group(diag, (/ diag%axesCu1%handles(1), diag%axesCu1%handles(2), z_axis /),    CS%axesCuz)
    call define_axes_group(diag, (/ diag%axesCv1%handles(1), diag%axesCv1%handles(2), z_axis /),    CS%axesCvz)
    call define_axes_group(diag, (/ diag%axesB1%handles(1),  diag%axesB1%handles(2),  zint_axis /), CS%axesBzi)
    call define_axes_group(diag, (/ diag%axesT1%handles(1),  diag%axesT1%handles(2),  zint_axis /), CS%axesTzi)
    call define_axes_group(diag, (/ diag%axesCu1%handles(1), diag%axesCu1%handles(2), zint_axis /), CS%axesCuzi)
    call define_axes_group(diag, (/ diag%axesCv1%handles(1), diag%axesCv1%handles(2), zint_axis /), CS%axesCvzi)
    call define_axes_group(diag, (/ z_axis /),    CS%axesZ)

    CS%id_u_z = register_diag_field('ocean_model_zold', 'u', CS%axesCuz, Time,    &
        'Zonal Velocity in Depth Space', 'm s-1',                     &
        missing_value=CS%missing_vel, cmor_field_name='uo', cmor_units='m s-1',&
        cmor_standard_name='sea_water_x_velocity', cmor_long_name='Sea Water X Velocity')
    if (CS%id_u_z>0) call safe_alloc_ptr(CS%u_z,IsdB,IedB,jsd,jed,CS%nk_zspace)

    CS%id_v_z = register_diag_field('ocean_model_zold', 'v', CS%axesCvz, Time,    &
        'Meridional Velocity in Depth Space', 'm s-1',                &
        missing_value=CS%missing_vel, cmor_field_name='vo', cmor_units='m s-1',&
        cmor_standard_name='sea_water_y_velocity', cmor_long_name='Sea Water Y Velocity')
    if (CS%id_v_z>0) call safe_alloc_ptr(CS%v_z,isd,ied,JsdB,JedB,CS%nk_zspace)

    CS%id_uh_z = register_diag_field('ocean_model_zold', 'uh', CS%axesCuz, Time,    &
        'Zonal Mass Transport (including SGS param) in Depth Space', flux_units, &
        missing_value=CS%missing_trans)
    if (CS%id_uh_z>0) call safe_alloc_ptr(CS%uh_z,IsdB,IedB,jsd,jed,CS%nk_zspace)

    CS%id_vh_z = register_diag_field('ocean_model_zold', 'vh', CS%axesCvz, Time,        &
        'Meridional Mass Transport (including SGS param) in Depth Space', flux_units,&
        missing_value=CS%missing_trans)
    if (CS%id_vh_z>0) call safe_alloc_ptr(CS%vh_z,isd,ied,JsdB,JedB,CS%nk_zspace)

  endif

end subroutine MOM_diag_to_Z_init

!> This subroutine reads the depths of the interfaces bounding the intended
!! layers from a NetCDF file.  If no appropriate file is found, -1 is returned
!! as the number of layers in the output file.  Also, a diag_manager axis is set
!! up with the same information as this axis.
subroutine get_Z_depths(depth_file, int_depth_name, int_depth, cell_depth_name, &
                        z_axis_index, edge_index, nk_out)
  character(len=*),   intent(in)  :: depth_file  !< The file to read for the depths
  character(len=*),   intent(in)  :: int_depth_name !< The interface depth variable name
  real, dimension(:), pointer     :: int_depth   !< A pointer that will be allocated and
                                                 !! returned with the interface depths
  character(len=*),   intent(in)  :: cell_depth_name !< The cell-center depth variable name
  integer,            intent(out) :: z_axis_index !< The cell-center z-axis diagnostic index handle
  integer,            intent(out) :: edge_index  !< The interface z-axis diagnostic index handle
  integer,            intent(out) :: nk_out      !< The number of layers in the output grid
  ! Local variables
  real, allocatable   :: cell_depth(:)
  character (len=200) :: units, long_name
  integer :: ncid, status, intid, intvid, layid, layvid, k, ni

  nk_out = -1

  status = NF90_OPEN(depth_file, NF90_NOWRITE, ncid)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"MOM_diag_to_Z get_Z_depths: "//&
        " Difficulties opening "//trim(depth_file)//" - "//&
        trim(NF90_STRERROR(status)))
    nk_out = -1 ; return
  endif

  status = NF90_INQ_DIMID(ncid, int_depth_name, intid)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"MOM_diag_to_Z get_Z_depths: "//&
      trim(NF90_STRERROR(status))//" Getting ID of dimension "//&
      trim(int_depth_name)//" in "//trim(depth_file))
    nk_out = -1 ; return
  endif

  status = nf90_Inquire_Dimension(ncid, intid, len=ni)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"MOM_diag_to_Z get_Z_depths: "//&
      trim(NF90_STRERROR(status))//" Getting number of interfaces of "//&
      trim(int_depth_name)//" in "//trim(depth_file))
    nk_out = -1 ; return
  endif

  if (ni < 2) then
    call MOM_error(WARNING,"MOM_diag_to_Z get_Z_depths: "//&
        "At least two interface depths must be specified in "//trim(depth_file))
    nk_out = -1 ; return
  endif

  status = NF90_INQ_DIMID(ncid, cell_depth_name, layid)
  if (status /= NF90_NOERR) call MOM_error(WARNING,"MOM_diag_to_Z get_Z_depths: "//&
      trim(NF90_STRERROR(status))//" Getting ID of dimension "//&
      trim(cell_depth_name)//" in "//trim(depth_file))

  status = nf90_Inquire_Dimension(ncid, layid, len=nk_out)
  if (status /= NF90_NOERR) call MOM_error(WARNING,"MOM_diag_to_Z get_Z_depths: "//&
      trim(NF90_STRERROR(status))//" Getting number of interfaces of "//&
      trim(cell_depth_name)//" in "//trim(depth_file))

  if (ni /= nk_out+1) then
    call MOM_error(WARNING,"MOM_diag_to_Z get_Z_depths: "//&
        "The interface depths must have one more point than cell centers in "//&
        trim(depth_file))
    nk_out = -1 ; return
  endif

  allocate(int_depth(nk_out+1))
  allocate(cell_depth(nk_out))

  status = NF90_INQ_VARID(ncid, int_depth_name, intvid)
  if (status /= NF90_NOERR) call MOM_error(FATAL,"MOM_diag_to_Z get_Z_depths: "//&
      trim(NF90_STRERROR(status))//" Getting ID of variable "//&
      trim(int_depth_name)//" in "//trim(depth_file))
  status = NF90_GET_VAR(ncid, intvid, int_depth)
  if (status /= NF90_NOERR) call MOM_error(FATAL,"MOM_diag_to_Z get_Z_depths: "//&
      trim(NF90_STRERROR(status))//" Reading variable "//&
      trim(int_depth_name)//" in "//trim(depth_file))
  status = NF90_GET_ATT(ncid, intvid, "units", units)
  if (status /= NF90_NOERR) units = "m"
  status = NF90_GET_ATT(ncid, intvid, "long_name", long_name)
  if (status /= NF90_NOERR) long_name = "Depth of edges"
  edge_index = diag_axis_init(int_depth_name, int_depth, units, 'z', &
                              long_name, direction=-1)

! Create an fms axis with the same data as the cell_depth array in the input file.
  status = NF90_INQ_VARID(ncid, cell_depth_name, layvid)
  if (status /= NF90_NOERR) call MOM_error(FATAL,"MOM_diag_to_Z get_Z_depths: "//&
      trim(NF90_STRERROR(status))//" Getting ID of variable "//&
      trim(cell_depth_name)//" in "//trim(depth_file))
  status = NF90_GET_VAR(ncid, layvid, cell_depth)
  if (status /= NF90_NOERR) call MOM_error(FATAL,"MOM_diag_to_Z get_Z_depths: "//&
      trim(NF90_STRERROR(status))//" Reading variable "//&
      trim(cell_depth_name)//" in "//trim(depth_file))
  status = NF90_GET_ATT(ncid, layvid, "units", units)
  if (status /= NF90_NOERR) units = "m"
  status = NF90_GET_ATT(ncid, layvid, "long_name", long_name)
  if (status /= NF90_NOERR) long_name = "Depth of cell center"

  z_axis_index = diag_axis_init(cell_depth_name, cell_depth, units, 'z',&
                                long_name, edges = edge_index, direction=-1)

  deallocate(cell_depth)

  status = nf90_close(ncid)

  ! Check the sign convention and change to the MOM "height" convention.
  if (int_depth(1) < int_depth(2)) then
    do k=1,nk_out+1 ; int_depth(k) = -1*int_depth(k) ; enddo
  endif

  ! Check for inversions in grid.
  do k=1,nk_out ; if (int_depth(k) < int_depth(k+1)) then
    call MOM_error(WARNING,"MOM_diag_to_Z get_Z_depths: "//&
        "Inverted interface depths in output grid in "//depth_file)
    nk_out = -1 ; deallocate(int_depth) ; return
  endif ; enddo

end subroutine get_Z_depths

!> Deallocate memory associated with the MOM_diag_to_Z module
subroutine MOM_diag_to_Z_end(CS)
  type(diag_to_Z_CS), pointer :: CS !< Control structure returned by a previous call to diag_to_Z_init.
  integer :: m

  if (associated(CS%u_z))   deallocate(CS%u_z)
  if (associated(CS%v_z))   deallocate(CS%v_z)
  if (associated(CS%Z_int)) deallocate(CS%Z_int)
  do m=1,CS%num_tr_used ;   deallocate(CS%tr_z(m)%p) ; enddo

  deallocate(CS)

end subroutine MOM_diag_to_Z_end

!> This subroutine registers a tracer to be output in depth space.
function ocean_register_diag_with_z(tr_ptr, vardesc_tr, G, Time, CS)
  type(ocean_grid_type), intent(in) :: G      !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                 target, intent(in) :: tr_ptr !< Tracer for translation to Z-space.
  type(vardesc),         intent(in) :: vardesc_tr !< Variable descriptor.
  type(time_type),       intent(in) :: Time   !< Current model time.
  type(diag_to_Z_CS),    pointer    :: CS     !< Control struct returned by a previous
                                              !! call to diag_to_Z_init.
  integer                           :: ocean_register_diag_with_z !< The retuned Z-space diagnostic ID
  ! Local variables
  type(vardesc) :: vardesc_z
  character(len=64) :: var_name         ! A variable's name.
  integer :: isd, ied, jsd, jed, nk, m, id_test
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nk = G%ke
  if (.not.associated(CS)) call MOM_error(FATAL, &
         "register_Z_tracer: Module must be initialized before it is used.")
  if (CS%nk_zspace<1) return

  if (CS%num_tr_used >= MAX_FIELDS_) then
    call MOM_error(WARNING,"ocean_register_diag_with_z:  Attempted to register and use "//&
                   "more than MAX_FIELDS_ z-space tracers via ocean_register_diag_with_z.")
    return
  endif

  ! register the layer tracer
  ocean_register_diag_with_z = ocean_register_diag(vardesc_tr, G, CS%diag, Time)

  ! copy layer tracer variable descriptor to a z-tracer descriptor
  ! change the name and layer information.
  vardesc_z = vardesc_tr
  call modify_vardesc(vardesc_z, z_grid="z", caller="ocean_register_diag_with_z")
  m = CS%num_tr_used + 1
  CS%missing_tr(m) = CS%missing_value ! This could be changed later, if desired.
  CS%id_tr(m) =  register_Z_diag(vardesc_z, CS, Time, CS%missing_tr(m))

  if (CS%nk_zspace > NO_ZSPACE) then
! There is a depth-space target file.
    if (CS%id_tr(m)>0) then
! Only allocate the tr_z field id there is a diag_table entry looking
! for it.
      CS%num_tr_used = m
      call safe_alloc_ptr(CS%tr_z(m)%p,isd,ied,jsd,jed,CS%nk_zspace)
!Can we do the following at this point?
! tr_ptr might not be allocated yet
      CS%tr_model(m)%p => tr_ptr
    endif
  else
! There is no depth-space target file but warn if a diag_table entry is
! present.
    call query_vardesc(vardesc_z, name=var_name, caller="ocean_register_diag_with_z")
    if (CS%id_tr(m)>0) call MOM_error(WARNING, &
        "ocean_register_diag_with_z: "//trim(var_name)// &
        " cannot be output without an appropriate depth-space target file.")
  endif

end function ocean_register_diag_with_z

!> Register a diagnostic to be output in depth space.
function register_Z_diag(var_desc, CS, day, missing)
  integer                        :: register_Z_diag !< The returned z-layer diagnostic index
  type(vardesc),      intent(in) :: var_desc !< A type with metadata for this diagnostic
  type(diag_to_Z_CS), pointer    :: CS   !< Control structure returned by
                                         !! previous call to diag_to_Z_init.
  type(time_type),    intent(in) :: day  !< The current model time
  real,               intent(in) :: missing !< The missing value for this diagnostic
  ! Local variables
  character(len=64) :: var_name         ! A variable's name.
  character(len=48) :: units            ! A variable's units.
  character(len=240) :: longname        ! A variable's longname.
  character(len=8) :: hor_grid, z_grid  ! Variable grid info.
  type(axes_grp), pointer :: axes => NULL()

  call query_vardesc(var_desc, name=var_name, units=units, longname=longname, &
                     hor_grid=hor_grid, z_grid=z_grid, caller="register_Zint_diag")

  ! Use the hor_grid and z_grid components of vardesc to determine the
  ! desired axes to register the diagnostic field for.
  select case (z_grid)

    case ("z")
      select case (hor_grid)
        case ("q")
          axes => CS%axesBz
        case ("h")
          axes => CS%axesTz
        case ("u")
          axes => CS%axesCuz
        case ("v")
          axes => CS%axesCvz
        case ("Bu")
          axes => CS%axesBz
        case ("T")
          axes => CS%axesTz
        case ("Cu")
          axes => CS%axesCuz
        case ("Cv")
          axes => CS%axesCvz
        case default
          call MOM_error(FATAL,&
            "register_Z_diag: unknown hor_grid component "//trim(hor_grid))
      end select

    case default
        call MOM_error(FATAL,&
          "register_Z_diag: unknown z_grid component "//trim(z_grid))
  end select

  register_Z_diag = register_diag_field("ocean_model_zold", trim(var_name), axes, &
        day, trim(longname), trim(units), missing_value=missing)

end function register_Z_diag

!> Register a diagnostic to be output at depth space interfaces
function register_Zint_diag(var_desc, CS, day, conversion)
  integer                        :: register_Zint_diag !< The returned z-interface diagnostic index
  type(vardesc),      intent(in) :: var_desc !< A type with metadata for this diagnostic
  type(diag_to_Z_CS), pointer    :: CS   !< Control structure returned by
                                         !! previous call to diag_to_Z_init.
  type(time_type),    intent(in) :: day  !< The current model time
  real,     optional, intent(in) :: conversion !< A value to multiply data by before writing to file
  ! Local variables
  character(len=64) :: var_name         ! A variable's name.
  character(len=48) :: units            ! A variable's units.
  character(len=240) :: longname        ! A variable's longname.
  character(len=8) :: hor_grid          ! Variable grid info.
  type(axes_grp), pointer :: axes => NULL()

  call query_vardesc(var_desc, name=var_name, units=units, longname=longname, &
                     hor_grid=hor_grid, caller="register_Zint_diag")

  if (CS%nk_zspace < 0) then
    register_Zint_diag = -1 ; return
  endif

  ! Use the hor_grid and z_grid components of vardesc to determine the
  ! desired axes to register the diagnostic field for.
  select case (hor_grid)
    case ("h")
      axes => CS%axesTzi
    case ("q")
      axes => CS%axesBzi
    case ("u")
      axes => CS%axesCuzi
    case ("v")
      axes => CS%axesCvzi
    case default
      call MOM_error(FATAL,&
        "register_Z_diag: unknown hor_grid component "//trim(hor_grid))
  end select

  register_Zint_diag = register_diag_field("ocean_model_zold", trim(var_name), &
        axes, day, trim(longname), trim(units), missing_value=CS%missing_value, &
        conversion=conversion)

end function register_Zint_diag

end module MOM_diag_to_Z
