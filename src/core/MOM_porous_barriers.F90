!> Module for calculating curve fit for porous topography.
!written by sjd
module MOM_porous_barriers

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,         only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_MODULE
use MOM_error_handler,     only : MOM_error, FATAL
use MOM_grid,              only : ocean_grid_type
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : thermo_var_ptrs, porous_barrier_type
use MOM_verticalGrid,      only : verticalGrid_type
use MOM_interface_heights, only : find_eta
use MOM_time_manager,      only : time_type
use MOM_diag_mediator,     only : register_diag_field, diag_ctrl, post_data
use MOM_file_parser,       only : param_file_type, get_param, log_version
use MOM_unit_scaling,      only : unit_scale_type
use MOM_debugging,         only : hchksum, uvchksum

implicit none ; private

public porous_widths_layer, porous_widths_interface, porous_barriers_init

#include <MOM_memory.h>

!> The control structure for the MOM_porous_barriers module
type, public :: porous_barrier_CS; private
  logical :: initialized = .false.  !< True if this control structure has been initialized.
  type(diag_ctrl), pointer :: &
      diag => Null()                !< A structure to regulate diagnostic output timing
  logical :: debug                  !< If true, write verbose checksums for debugging purposes.
  real    :: mask_depth             !< The depth shallower than which porous barrier is not applied [Z ~> m]
  integer :: eta_interp             !< An integer indicating how the interface heights at the velocity
                                    !! points are calculated. Valid values are given by the parameters
                                    !! defined below: MAX, MIN, ARITHMETIC and HARMONIC.
  integer :: answer_date            !< The vintage of the porous barrier weight function calculations.
                                    !! Values below 20220806 recover the old answers in which the layer
                                    !! averaged weights are not strictly limited by an upper-bound of 1.0 .
  !>@{ Diagnostic IDs
  integer :: id_por_layer_widthU = -1, id_por_layer_widthV = -1, &
             id_por_face_areaU = -1, id_por_face_areaV = -1
  !>@}
end type porous_barrier_CS

integer :: id_clock_porous_barrier !< CPU clock for porous barrier

!>@{ Enumeration values for eta interpolation schemes
integer, parameter :: ETA_INTERP_MAX   = 1
integer, parameter :: ETA_INTERP_MIN   = 2
integer, parameter :: ETA_INTERP_ARITH = 3
integer, parameter :: ETA_INTERP_HARM  = 4
character(len=20), parameter :: ETA_INTERP_MAX_STRING = "MAX"
character(len=20), parameter :: ETA_INTERP_MIN_STRING = "MIN"
character(len=20), parameter :: ETA_INTERP_ARITH_STRING = "ARITHMETIC"
character(len=20), parameter :: ETA_INTERP_HARM_STRING = "HARMONIC"
!>@}

contains

!> subroutine to assign porous barrier widths averaged over a layer
subroutine porous_widths_layer(h, tv, G, GV, US, pbv, CS, eta_bt)
  ! Note: eta_bt is not currently used
  type(ocean_grid_type),                      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in) :: tv  !< A structure pointing to various
                                                                !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: eta_bt !< optional barotropic variable
                                                                   !! used to dilate the layer thicknesses
                                                                   !! [H ~> m or kg m-2].
  type(porous_barrier_type),                  intent(inout) :: pbv !< porous barrier fractional cell metrics
  type(porous_barrier_CS),                    intent(in) :: CS     !< Control structure for porous barrier

  !local variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: eta_u ! Layer interface heights at u points [Z ~> m]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: eta_v ! Layer interface heights at v points [Z ~> m]
  real, dimension(SZIB_(G),SZJB_(G)) :: A_layer_prev ! Integral of fractional open width from the bottom
                                                     ! to the previous layer at u or v points [Z ~> m]
  logical, dimension(SZIB_(G),SZJB_(G)) :: do_I ! Booleans for calculation at u or v points
                                                ! updated while moving up layers
  real :: A_layer ! Integral of fractional open width from bottom to current layer [Z ~> m]
  real :: dz_min  ! The minimum layer thickness [Z ~> m]
  real :: dmask ! The depth below which porous barrier is not applied [Z ~> m]
  integer :: i, j, k, nk, is, ie, js, je, Isq, Ieq, Jsq, Jeq

  if (.not.CS%initialized) call MOM_error(FATAL, &
      "MOM_Porous_barrier: Module must be initialized before it is used.")

  call cpu_clock_begin(id_clock_porous_barrier)

  is = G%isc; ie = G%iec; js = G%jsc; je = G%jec; nk = GV%ke
  Isq = G%IscB; Ieq = G%IecB; Jsq = G%JscB; Jeq = G%JecB

  if (CS%answer_date < 20220806) then
    dmask = 0.0
  else
    dmask = CS%mask_depth
  endif

  call calc_eta_at_uv(eta_u, eta_v, CS%eta_interp, dmask, h, tv, G, GV, US)

  dz_min = GV%Angstrom_Z

  ! u-points
  do j=js,je ; do I=Isq,Ieq ; do_I(I,j) = .False. ; enddo ; enddo

  do j=js,je ; do I=Isq,Ieq ; if (G%porous_DavgU(I,j) < dmask) then
    call calc_por_layer(G%porous_DminU(I,j), G%porous_DmaxU(I,j), G%porous_DavgU(I,j), &
                        eta_u(I,j,nk+1), A_layer_prev(I,j), do_I(I,j))
  endif ; enddo ; enddo

  if (CS%answer_date < 20220806) then
    do k=nk,1,-1 ; do j=js,je ; do I=Isq,Ieq ; if (G%porous_DavgU(I,j) < dmask) then
      call calc_por_layer(G%porous_DminU(I,j), G%porous_DmaxU(I,j), G%porous_DavgU(I,j), &
                          eta_u(I,j,K), A_layer, do_I(I,j))
      if (eta_u(I,j,K) - eta_u(I,j,K+1) > 0.0) then
        pbv%por_face_areaU(I,j,k) = (A_layer - A_layer_prev(I,j)) / (eta_u(I,j,K) - eta_u(I,j,K+1))
      else
        pbv%por_face_areaU(I,j,k) = 0.0
      endif
      A_layer_prev(I,j) = A_layer
    endif ; enddo ; enddo ; enddo
  else
    do k=nk,1,-1 ; do j=js,je ; do I=Isq,Ieq
      if (do_I(I,j)) then
        call calc_por_layer(G%porous_DminU(I,j), G%porous_DmaxU(I,j), G%porous_DavgU(I,j), &
                            eta_u(I,j,K), A_layer, do_I(I,j))
        if (eta_u(I,j,K) - (eta_u(I,j,K+1)+dz_min) > 0.0) then
          pbv%por_face_areaU(I,j,k) = min(1.0, (A_layer - A_layer_prev(I,j)) / (eta_u(I,j,K) - eta_u(I,j,K+1)))
        else
          pbv%por_face_areaU(I,j,k) = 0.0 ! use calc_por_interface() might be a better choice
        endif
        A_layer_prev(I,j) = A_layer
      else
        pbv%por_face_areaU(I,j,k) = 1.0
      endif
    enddo ; enddo ; enddo
  endif

  ! v-points
  do J=Jsq,Jeq ; do i=is,ie; do_I(i,J) = .False. ; enddo ; enddo

  do J=Jsq,Jeq ; do i=is,ie ; if (G%porous_DavgV(i,J) < dmask) then
    call calc_por_layer(G%porous_DminV(i,J), G%porous_DmaxV(i,J), G%porous_DavgV(i,J), &
                        eta_v(i,J,nk+1), A_layer_prev(i,J), do_I(i,J))
  endif ; enddo ; enddo

  if (CS%answer_date < 20220806) then
    do k=nk,1,-1 ; do J=Jsq,Jeq ; do i=is,ie ; if (G%porous_DavgV(i,J) < dmask) then
      call calc_por_layer(G%porous_DminV(i,J), G%porous_DmaxV(i,J), G%porous_DavgV(i,J), &
                          eta_v(i,J,K), A_layer, do_I(i,J))
      if (eta_v(i,J,K) - eta_v(i,J,K+1) > 0.0) then
        pbv%por_face_areaV(i,J,k) = (A_layer - A_layer_prev(i,J)) / (eta_v(i,J,K) - eta_v(i,J,K+1))
      else
        pbv%por_face_areaV(i,J,k) = 0.0
      endif
      A_layer_prev(i,J) = A_layer
    endif ; enddo ; enddo ; enddo
  else
    do k=nk,1,-1 ; do J=Jsq,Jeq ; do i=is,ie
      if (do_I(i,J)) then
        call calc_por_layer(G%porous_DminV(i,J), G%porous_DmaxV(i,J), G%porous_DavgV(i,J), &
                            eta_v(i,J,K), A_layer, do_I(i,J))
        if (eta_v(i,J,K) - (eta_v(i,J,K+1)+dz_min) > 0.0) then
          pbv%por_face_areaV(i,J,k) = min(1.0, (A_layer - A_layer_prev(i,J)) / (eta_v(i,J,K) - eta_v(i,J,K+1)))
        else
          pbv%por_face_areaV(i,J,k) = 0.0 ! use calc_por_interface() might be a better choice
        endif
        A_layer_prev(i,J) = A_layer
      else
        pbv%por_face_areaV(i,J,k) = 1.0
      endif
    enddo ; enddo ; enddo
  endif

  if (CS%debug) then
    call uvchksum("Interface height used by porous barrier for layer weights", &
                  eta_u, eta_v, G%HI, haloshift=0, scalar_pair=.true.)
    call uvchksum("Porous barrier layer-averaged weights: por_face_area[UV]", &
                  pbv%por_face_areaU, pbv%por_face_areaV, G%HI, haloshift=0, &
                  scalar_pair=.true.)
  endif

  if (CS%id_por_face_areaU > 0) call post_data(CS%id_por_face_areaU, pbv%por_face_areaU, CS%diag)
  if (CS%id_por_face_areaV > 0) call post_data(CS%id_por_face_areaV, pbv%por_face_areaV, CS%diag)

  call cpu_clock_end(id_clock_porous_barrier)
end subroutine porous_widths_layer

!> subroutine to assign porous barrier widths at the layer interfaces
subroutine porous_widths_interface(h, tv, G, GV, US, pbv, CS, eta_bt)
  ! Note: eta_bt is not currently used
  type(ocean_grid_type),                      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in) :: tv  !< A structure pointing to various
                                                                !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: eta_bt !< optional barotropic variable
                                                                   !! used to dilate the layer thicknesses
                                                                   !! [H ~> m or kg m-2].
  type(porous_barrier_type),                  intent(inout) :: pbv  !< porous barrier fractional cell metrics
  type(porous_barrier_CS),                    intent(in) :: CS !< Control structure for porous barrier

  !local variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: eta_u ! Layer interface height at u points [Z ~> m]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: eta_v ! Layer interface height at v points [Z ~> m]
  logical, dimension(SZIB_(G),SZJB_(G)) :: do_I ! Booleans for calculation at u or v points
                                                ! updated while moving up layers
  real :: dmask ! The depth below which porous barrier is not applied [Z ~> m]
  integer :: i, j, k, nk, is, ie, js, je, Isq, Ieq, Jsq, Jeq

  if (.not.CS%initialized) call MOM_error(FATAL, &
      "MOM_Porous_barrier: Module must be initialized before it is used.")

  call cpu_clock_begin(id_clock_porous_barrier)

  is = G%isc; ie = G%iec; js = G%jsc; je = G%jec; nk = GV%ke
  Isq = G%IscB; Ieq = G%IecB; Jsq = G%JscB; Jeq = G%JecB

  if (CS%answer_date < 20220806) then
    dmask = 0.0
  else
    dmask = CS%mask_depth
  endif

  call calc_eta_at_uv(eta_u, eta_v, CS%eta_interp, dmask, h, tv, G, GV, US)

  ! u-points
  do j=js,je ; do I=Isq,Ieq
    do_I(I,j) = .False.
    if (G%porous_DavgU(I,j) < dmask) do_I(I,j) = .True.
  enddo ; enddo

  if (CS%answer_date < 20220806) then
    do K=1,nk+1 ; do j=js,je ; do I=Isq,Ieq ; if (G%porous_DavgU(I,j) < dmask) then
      call calc_por_interface(G%porous_DminU(I,j), G%porous_DmaxU(I,j), G%porous_DavgU(I,j), &
                              eta_u(I,j,K), pbv%por_layer_widthU(I,j,K), do_I(I,j))
   endif ; enddo ; enddo ; enddo
  else
    do K=1,nk+1 ; do j=js,je ; do I=Isq,Ieq
      if (do_I(I,j)) then
        call calc_por_interface(G%porous_DminU(I,j), G%porous_DmaxU(I,j), G%porous_DavgU(I,j), &
                                eta_u(I,j,K), pbv%por_layer_widthU(I,j,K), do_I(I,j))
      else
        pbv%por_layer_widthU(I,j,K) = 1.0
      endif
    enddo ; enddo ; enddo
  endif

  ! v-points
  do J=Jsq,Jeq ; do i=is,ie
    do_I(i,J) = .False.
    if (G%porous_DavgV(i,J) < dmask) do_I(i,J) = .True.
  enddo ; enddo

  if (CS%answer_date < 20220806) then
    do K=1,nk+1 ; do J=Jsq,Jeq ; do i=is,ie ; if (G%porous_DavgV(i,J) < dmask) then
      call calc_por_interface(G%porous_DminV(i,J), G%porous_DmaxV(i,J), G%porous_DavgV(i,J), &
                              eta_v(i,J,K), pbv%por_layer_widthV(i,J,K), do_I(i,J))
    endif ; enddo ; enddo ; enddo
  else
    do K=1,nk+1 ; do J=Jsq,Jeq ; do i=is,ie
      if (do_I(i,J)) then
        call calc_por_interface(G%porous_DminV(i,J), G%porous_DmaxV(i,J), G%porous_DavgV(i,J), &
                                eta_v(i,J,K), pbv%por_layer_widthV(i,J,K), do_I(i,J))
      else
        pbv%por_layer_widthV(i,J,K) = 1.0
      endif
    enddo ; enddo ; enddo
  endif

  if (CS%debug) then
    call uvchksum("Interface height used by porous barrier for interface weights", &
                  eta_u, eta_v, G%HI, haloshift=0, scalar_pair=.true.)
    call uvchksum("Porous barrier weights at the layer-interface: por_layer_width[UV]", &
                  pbv%por_layer_widthU, pbv%por_layer_widthV, G%HI, &
                  haloshift=0, scalar_pair=.true.)
  endif

  if (CS%id_por_layer_widthU > 0) call post_data(CS%id_por_layer_widthU, pbv%por_layer_widthU, CS%diag)
  if (CS%id_por_layer_widthV > 0) call post_data(CS%id_por_layer_widthV, pbv%por_layer_widthV, CS%diag)

  call cpu_clock_end(id_clock_porous_barrier)
end subroutine porous_widths_interface

subroutine calc_eta_at_uv(eta_u, eta_v, interp, dmask, h, tv, G, GV, US, eta_bt)
  !variables needed to call find_eta
  type(ocean_grid_type),                        intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                      intent(in) :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),                        intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                        intent(in) :: tv  !< A structure pointing to various
                                                                  !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G)), optional,   intent(in) :: eta_bt !< optional barotropic variable
                                                                   !! used to dilate the layer thicknesses
                                                                   !! [H ~> m or kg m-2].
  real,                                         intent(in) :: dmask !< The depth shallower than which
                                                                    !! porous barrier is not applied [Z ~> m]
  integer,                                      intent(in) :: interp !< eta interpolation method
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(out) :: eta_u !< Layer interface heights at u points [Z ~> m]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(out) :: eta_v !< Layer interface heights at v points [Z ~> m]

  ! local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: eta ! Layer interface heights [Z ~> m].
  real :: dz_neglect ! A negligible height difference [Z ~> m]
  integer :: i, j, k, nk, is, ie, js, je, Isq, Ieq, Jsq, Jeq

  is = G%isc; ie = G%iec; js = G%jsc; je = G%jec; nk = GV%ke
  Isq = G%IscB; Ieq = G%IecB; Jsq = G%JscB; Jeq = G%JecB

  ! currently no treatment for using optional find_eta arguments if present
  call find_eta(h, tv, G, GV, US, eta, halo_size=1)

  dz_neglect = GV%dZ_subroundoff

  do K=1,nk+1
    do j=js,je ; do I=Isq,Ieq ; eta_u(I,j,K) = dmask ; enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie ; eta_v(i,J,K) = dmask ; enddo ; enddo
  enddo

  select case (interp)
    case (ETA_INTERP_MAX)   ! The shallower interface height
      do K=1,nk+1
        do j=js,je ; do I=Isq,Ieq ; if (G%porous_DavgU(I,j) < dmask) then
          eta_u(I,j,K) = max(eta(i,j,K), eta(i+1,j,K))
        endif ; enddo ; enddo
        do J=Jsq,Jeq ; do i=is,ie ; if (G%porous_DavgV(i,J) < dmask) then
          eta_v(i,J,K) = max(eta(i,j,K), eta(i,j+1,K))
        endif ; enddo ; enddo
      enddo
    case (ETA_INTERP_MIN)   ! The deeper interface height
      do K=1,nk+1
        do j=js,je ; do I=Isq,Ieq ; if (G%porous_DavgU(I,j) < dmask) then
          eta_u(I,j,K) = min(eta(i,j,K), eta(i+1,j,K))
        endif ; enddo ; enddo
        do J=Jsq,Jeq ; do i=is,ie ; if (G%porous_DavgV(i,J) < dmask) then
          eta_v(i,J,K) = min(eta(i,j,K), eta(i,j+1,K))
        endif ; enddo ; enddo
      enddo
    case (ETA_INTERP_ARITH) ! Arithmetic mean
      do K=1,nk+1
        do j=js,je ; do I=Isq,Ieq ; if (G%porous_DavgU(I,j) < dmask) then
          eta_u(I,j,K) = 0.5 * (eta(i,j,K) + eta(i+1,j,K))
        endif ; enddo ; enddo
        do J=Jsq,Jeq ; do i=is,ie ; if (G%porous_DavgV(i,J) < dmask) then
          eta_v(i,J,K) = 0.5 * (eta(i,j,K) + eta(i,j+1,K))
        endif ; enddo ; enddo
      enddo
    case (ETA_INTERP_HARM)  ! Harmonic mean
      do K=1,nk+1
        do j=js,je ; do I=Isq,Ieq ; if (G%porous_DavgU(I,j) < dmask) then
          eta_u(I,j,K) = 2.0 * (eta(i,j,K) * eta(i+1,j,K)) / (eta(i,j,K) + eta(i+1,j,K) + dz_neglect)
        endif ; enddo ; enddo
        do J=Jsq,Jeq ; do i=is,ie ; if (G%porous_DavgV(i,J) < dmask) then
          eta_v(i,J,K) = 2.0 * (eta(i,j,K) * eta(i,j+1,K)) / (eta(i,j,K) + eta(i,j+1,K) + dz_neglect)
        endif ; enddo ; enddo
      enddo
    case default
      call MOM_error(FATAL, "porous_widths::calc_eta_at_uv: "//&
                     "invalid value for eta interpolation method.")
  end select
end subroutine calc_eta_at_uv

!> subroutine to calculate the profile fit (the three parameter fit from Adcroft 2013)
! of the open face area fraction below a certain depth (eta_layer) in a column
subroutine calc_por_layer(D_min, D_max, D_avg, eta_layer, A_layer, do_next)
  real,    intent(in)  :: D_min     !< minimum topographic height (deepest) [Z ~> m]
  real,    intent(in)  :: D_max     !< maximum topographic height (shallowest) [Z ~> m]
  real,    intent(in)  :: D_avg     !< mean topographic height [Z ~> m]
  real,    intent(in)  :: eta_layer !< height of interface [Z ~> m]
  real,    intent(out) :: A_layer   !< frac. open face area of below eta_layer [Z ~> m]
  logical, intent(out) :: do_next   !< False if eta_layer>D_max

  ! local variables
  real :: m      ! convenience constant for fit [nondim]
  real :: zeta   ! normalized vertical coordinate [nondim]

  do_next = .True.
  if (eta_layer <= D_min) then
    A_layer = 0.0
  elseif (eta_layer > D_max) then
    A_layer = eta_layer - D_avg
    do_next = .False.
  else
    m = (D_avg - D_min) / (D_max - D_min)
    zeta = (eta_layer - D_min) / (D_max - D_min)
    if (m < 0.5) then
      A_layer = (D_max - D_min) * ((1.0 - m) * zeta**(1.0 / (1.0 - m)))
    elseif (m == 0.5) then
      A_layer = (D_max - D_min) * (0.5 * zeta * zeta)
    else
      A_layer = (D_max - D_min) * (zeta - m + m * ((1.0 - zeta)**(1.0 / m)))
    endif
  endif
end subroutine calc_por_layer

!> subroutine to calculate the profile fit (the three parameter fit from Adcroft 2013)
! of the open interface fraction at a certain depth (eta_layer) in a column
subroutine calc_por_interface(D_min, D_max, D_avg, eta_layer, w_layer, do_next)
  real,    intent(in)  :: D_min     !< minimum topographic height (deepest) [Z ~> m]
  real,    intent(in)  :: D_max     !< maximum topographic height (shallowest) [Z ~> m]
  real,    intent(in)  :: D_avg     !< mean topographic height [Z ~> m]
  real,    intent(in)  :: eta_layer !< height of interface [Z ~> m]
  real,    intent(out) :: w_layer   !< frac. open interface width at eta_layer [nondim]
  logical, intent(out) :: do_next   !< False if eta_layer>D_max

  ! local variables
  real :: m, a     ! convenience constants for fit [nondim]
  real :: zeta     ! normalized vertical coordinate [nondim]

  do_next = .True.
  if (eta_layer <= D_min) then
    w_layer = 0.0
  elseif (eta_layer > D_max) then
    w_layer = 1.0
    do_next = .False.
  else  ! The following option could be refactored for stability and efficiency (with fewer divisions)
    m = (D_avg - D_min) / (D_max - D_min)
    a = (1.0 - m) / m
    zeta = (eta_layer - D_min) / (D_max - D_min)
    if (m < 0.5) then
      w_layer = zeta**(1.0 / a)
      ! Note that this would be safer and more efficent if it were rewritten as:
      ! w_layer = zeta**( (D_avg - D_min) / (D_max - D_avg) )
    elseif (m == 0.5) then
      w_layer = zeta
    else
      w_layer = 1.0 - (1.0 - zeta)**a
    endif
  endif
end subroutine calc_por_interface

subroutine porous_barriers_init(Time, GV, US, param_file, diag, CS)
  type(time_type),         intent(in)    :: Time       !< Current model time
  type(verticalGrid_type), intent(in)    :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< structure indicating parameter file to parse
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics control structure
  type(porous_barrier_CS), intent(inout) :: CS         !< Module control structure

  ! local variables
  character(len=40) :: mdl = "MOM_porous_barriers"  ! This module's name.
  character(len=20) :: interp_method ! String storing eta interpolation method
  integer :: default_answer_date ! Global answer date
  !> This include declares and sets the variable "version".
# include "version_variable.h"

  CS%initialized = .true.
  CS%diag => diag

  call log_version(param_file, mdl, version, "", log_to_all=.true., layout=.false., &
                   debugging=.false.)
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "PORBAR_ANSWER_DATE", CS%answer_date, &
                 "The vintage of the porous barrier weight function calculations.  Values below "//&
                 "20220806 recover the old answers in which the layer averaged weights are not "//&
                 "strictly limited by an upper-bound of 1.0 .", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) CS%answer_date = max(CS%answer_date, 20230701)
  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false.)
  call get_param(param_file, mdl, "PORBAR_MASKING_DEPTH", CS%mask_depth, &
                 "If the effective average depth at the velocity cell is shallower than this "//&
                 "number, then porous barrier is not applied at that location.  "//&
                 "PORBAR_MASKING_DEPTH is assumed to be positive below the sea surface.", &
                 units="m", default=0.0, scale=US%m_to_Z)
  ! The sign needs to be inverted to be consistent with the sign convention of Davg_[UV]
  CS%mask_depth = -CS%mask_depth
  call get_param(param_file, mdl, "PORBAR_ETA_INTERP", interp_method, &
                 "A string describing the method that decides how the "//&
                 "interface heights at the velocity points are calculated. "//&
                 "Valid values are:\n"//&
                 "\t MAX (the default) - maximum of the adjacent cells \n"//&
                 "\t MIN - minimum of the adjacent cells \n"//&
                 "\t ARITHMETIC - arithmetic mean of the adjacent cells \n"//&
                 "\t HARMONIC - harmonic mean of the adjacent cells \n", &
                 default=ETA_INTERP_MAX_STRING)
  select case (interp_method)
    case (ETA_INTERP_MAX_STRING) ; CS%eta_interp = ETA_INTERP_MAX
    case (ETA_INTERP_MIN_STRING) ; CS%eta_interp = ETA_INTERP_MIN
    case (ETA_INTERP_ARITH_STRING) ; CS%eta_interp = ETA_INTERP_ARITH
    case (ETA_INTERP_HARM_STRING) ; CS%eta_interp = ETA_INTERP_HARM
    case default
      call MOM_error(FATAL, "porous_barriers_init: Unrecognized setting "// &
            "#define PORBAR_ETA_INTERP "//trim(interp_method)//" found in input file.")
  end select

  CS%id_por_layer_widthU = register_diag_field('ocean_model', 'por_layer_widthU', diag%axesCui, Time, &
     'Porous barrier open width fraction (at the layer interfaces) of the u-faces', 'nondim')
  CS%id_por_layer_widthV = register_diag_field('ocean_model', 'por_layer_widthV', diag%axesCvi, Time, &
     'Porous barrier open width fraction (at the layer interfaces) of the v-faces', 'nondim')
  CS%id_por_face_areaU = register_diag_field('ocean_model', 'por_face_areaU', diag%axesCuL, Time, &
     'Porous barrier open area fraction (layer averaged) of U-faces', 'nondim')
  CS%id_por_face_areaV = register_diag_field('ocean_model', 'por_face_areaV', diag%axesCvL, Time, &
     'Porous barrier open area fraction (layer averaged) of V-faces', 'nondim')

  id_clock_porous_barrier = cpu_clock_id('(Ocean porous barrier)', grain=CLOCK_MODULE)
end subroutine

end module MOM_porous_barriers
