!> Module for calculating curve fit for porous topography.
!written by sjd
module MOM_porous_barriers

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,         only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_MODULE
use MOM_error_handler,     only : MOM_error, FATAL
use MOM_grid,              only : ocean_grid_type
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : thermo_var_ptrs, porous_barrier_ptrs
use MOM_verticalGrid,      only : verticalGrid_type
use MOM_interface_heights, only : find_eta
use MOM_time_manager,      only : time_type
use MOM_diag_mediator,     only : register_diag_field, diag_ctrl, post_data
use MOM_file_parser,       only : param_file_type, get_param, log_version
use MOM_unit_scaling,      only : unit_scale_type
use MOM_debugging,         only : hchksum, uvchksum

implicit none ; private

public porous_widths, porous_barriers_init

#include <MOM_memory.h>

type, public :: porous_barrier_CS; private
  logical :: initialized = .false.  !< True if this control structure has been initialized.
  type(diag_ctrl), pointer :: diag => Null()   !< A structure to regulate diagnostic output timing
  logical :: debug  !< If true, write verbose checksums for debugging purposes.
  real :: mask_depth  !< The depth below which porous barrier is not applied.
  integer :: eta_interp !< An integer indicating how the interface heights at the velocity points
                        !! are calculated. Valid values are given by the parameters defined below:
                        !! MAX, MIN, ARITHMETIC and HARMONIC.
  integer :: id_por_layer_widthU = -1, id_por_layer_widthV = -1, id_por_face_areaU = -1, id_por_face_areaV = -1
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

!> subroutine to assign cell face areas and layer widths for porous topography
subroutine porous_widths(h, tv, G, GV, US, pbv, CS, eta_bt, halo_size, eta_to_m)
  !eta_bt, halo_size, eta_to_m not currently used
  !variables needed to call find_eta
  type(ocean_grid_type),                      intent(in)  :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)  :: tv     !< A structure pointing to various
                                                                    !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)  :: eta_bt !< optional barotropic
             !! variable that gives the "correct" free surface height (Boussinesq) or total water
             !! column mass per unit area (non-Boussinesq).  This is used to dilate the layer.
             !! thicknesses when calculating interfaceheights [H ~> m or kg m-2].
  integer,                          optional, intent(in)  :: halo_size !< width of halo points on
                                                                       !! which to calculate eta.

  real,                             optional, intent(in)  :: eta_to_m  !< The conversion factor from
             !! the units of eta to m; by default this is US%Z_to_m.
  type(porous_barrier_ptrs),           intent(inout) :: pbv  !< porous barrier fractional cell metrics
  type(porous_barrier_CS), intent(in) :: CS

  !local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1):: eta ! Layer interface heights [Z ~> m or 1/eta_to_m].
  real, dimension(SZIB_(G),SZJB_(G)) :: A_layer_prev ! Integral of fractional open width from the bottom
                                                     ! to the previous layer at u or v points [Z ~> m]
  real, dimension(SZIB_(G),SZJB_(G)) :: eta_prev ! Layter interface height of the previous layer
                                                 ! at u or v points [Z ~> m]
  real :: A_layer, & ! Integral of fractional open width from bottom to current layer [Z ~> m]
          eta_s ! Layer height used for fit [Z ~> m]
  real :: Z_to_eta, H_to_eta ! Unit conversion factors for eta.
  real :: h_neglect, & ! ! Negligible thicknesses, often [Z ~> m]
          h_min ! ! The minimum layer thickness, often [Z ~> m]
  real :: dmask ! The depth below which porous barrier is not applied.
  integer :: i, j, k, nk, is, ie, js, je, Isq, Ieq, Jsq, Jeq

  if (.not.CS%initialized) call MOM_error(FATAL, &
      "MOM_Porous_barrier: Module must be initialized before it is used.")

  call cpu_clock_begin(id_clock_porous_barrier)

  is = G%isc; ie = G%iec; js = G%jsc; je = G%jec; nk = GV%ke
  Isq = G%IscB; Ieq = G%IecB; Jsq = G%JscB; Jeq = G%JecB

  dmask = CS%mask_depth

  !eta is zero at surface and decreases downward
  !currently no treatment for using optional find_eta arguments if present
  call find_eta(h, tv, G, GV, US, eta, halo_size=1)

  Z_to_eta = 1.0 ; if (present(eta_to_m)) Z_to_eta = US%Z_to_m / eta_to_m
  H_to_eta = GV%H_to_m * US%m_to_Z * Z_to_eta
  h_neglect = GV%H_subroundoff * H_to_eta
  h_min = GV%Angstrom_H * H_to_eta

  ! u-points
  do j=js,je ; do I=Isq,Ieq ; if (G%porous_DavgU(I,j) < dmask) then
    eta_prev(I,j) = eta_at_uv(eta(i,j,nk+1), eta(i+1,j,nk+1), CS%eta_interp, h_neglect)
    ! eta_prev(I,j) = max(eta(i,j,nk+1), eta(i+1,j,nk+1))
    call calc_por_layer(G%porous_DminU(I,j), G%porous_DmaxU(I,j), G%porous_DavgU(I,j), &
                       eta_prev(I,j), pbv%por_layer_widthU(I,j,nk+1), A_layer_prev(I,j))
  endif ; enddo ; enddo

  do k = nk,1,-1; do j=js,je; do I=Isq,Ieq ; if (G%porous_DavgU(I,j) < dmask) then
    eta_s = eta_at_uv(eta(i,j,K), eta(i+1,j,K), CS%eta_interp, h_neglect)
    ! eta_s = max(eta(i,j,K), eta(i+1,j,K))
    if ((eta_s - eta_prev(I,j)) > 0.0) then
      call calc_por_layer(G%porous_DminU(I,j), G%porous_DmaxU(I,j), G%porous_DavgU(I,j), &
                          eta_s, pbv%por_layer_widthU(I,j,K), A_layer)
      pbv%por_face_areaU(I,j,k) = (A_layer - A_layer_prev(I,j)) / (eta_s - eta_prev(I,j))
    else
      pbv%por_face_areaU(I,j,k) = 0.0
    endif
    eta_prev(I,j) = eta_s
    A_layer_prev(I,j) = A_layer
  endif ; enddo ; enddo ; enddo

  ! v-points
  do J=Jsq,Jeq ; do i=is,ie ; if (G%porous_DavgV(i,J) < dmask) then
    eta_prev(i,J) = eta_at_uv(eta(i,j,nk+1), eta(i,j+1,nk+1), CS%eta_interp, h_neglect)
    ! eta_prev(i,J) = max(eta(i,j,nk+1), eta(i,j+1,nk+1))
    call calc_por_layer(G%porous_DminV(i,J), G%porous_DmaxV(i,J), G%porous_DavgV(i,J), &
                        eta_prev(i,J), pbv%por_layer_widthV(i,J,nk+1), A_layer_prev(i,J))
  endif ; enddo ; enddo

  do k = nk,1,-1; do J=Jsq,Jeq ; do i=is,ie ;if (G%porous_DavgV(i,J) < dmask) then
    eta_s = eta_at_uv(eta(i,j,K), eta(i,j+1,K), CS%eta_interp, h_neglect)
    ! eta_s = max(eta(i,j,K), eta(i,j+1,K))
    if ((eta_s - eta_prev(i,J)) > 0.0) then
      call calc_por_layer(G%porous_DminV(i,J), G%porous_DmaxV(i,J), G%porous_DavgV(i,J), &
                         eta_s, pbv%por_layer_widthV(i,J,K), A_layer)
      pbv%por_face_areaV(i,J,k) = (A_layer - A_layer_prev(i,J)) / (eta_s - eta_prev(i,J))
    else
      pbv%por_face_areaV(i,J,k) = 0.0
    endif
    eta_prev(i,J) = eta_s
    A_layer_prev(i,J) = A_layer
  endif ; enddo ; enddo ; enddo

  if (CS%debug) then
    call hchksum(eta, "Interface height used by porous barrier", G%HI, haloshift=0, scale=GV%H_to_m)
    call uvchksum("Porous barrier weights at the layer-interface: por_layer_width[UV]", &
                   pbv%por_layer_widthU, pbv%por_layer_widthV, G%HI, haloshift=0)
    call uvchksum("Porous barrier layer-averaged weights: por_face_area[UV]", &
                   pbv%por_face_areaU, pbv%por_face_areaV, G%HI, haloshift=0)
  endif

  if (CS%id_por_layer_widthU > 0) call post_data(CS%id_por_layer_widthU, pbv%por_layer_widthU, CS%diag)
  if (CS%id_por_layer_widthV > 0) call post_data(CS%id_por_layer_widthV, pbv%por_layer_widthV, CS%diag)
  if (CS%id_por_face_areaU > 0) call post_data(CS%id_por_face_areaU, pbv%por_face_areaU, CS%diag)
  if (CS%id_por_face_areaV > 0) call post_data(CS%id_por_face_areaV, pbv%por_face_areaV, CS%diag)

  call cpu_clock_end(id_clock_porous_barrier)
end subroutine porous_widths

!> subroutine to calculate the profile fit (the three parameter fit from Adcroft 2013)
! for a single layer in a column
subroutine calc_por_layer(D_min, D_max, D_avg, eta_layer, w_layer, A_layer)
  real, intent(in)  :: D_min !< minimum topographic height [Z ~> m]
  real, intent(in)  :: D_max !< maximum topographic height [Z ~> m]
  real, intent(in)  :: D_avg !< mean topographic height [Z ~> m]
  real, intent(in)  :: eta_layer !< height of interface [Z ~> m]
  real, intent(out) :: w_layer !< frac. open interface width of current layer [nondim]
  real, intent(out) :: A_layer !< frac. open face area of current layer [Z ~> m]

  ! local variables
  real :: m, a, &  ! convenience constant for fit [nondim]
          zeta, &  ! normalized vertical coordinate [nondim]
          psi, &   ! fractional width of layer between D_min and D_max [nondim]
          psi_int  ! integral of psi from 0 to zeta

  if (eta_layer <= D_min) then
    w_layer = 0.0
    A_layer = 0.0
  elseif (eta_layer > D_max) then
    w_layer = 1.0
    A_layer = eta_layer - D_avg
  else
    m = (D_avg - D_min) / (D_max - D_min)
    a = (1.0 - m) / m
    zeta = (eta_layer - D_min) / (D_max - D_min)
    if (m < 0.5) then
      psi = zeta**(1.0 / a)
      psi_int = (1.0 - m) * zeta**(1.0 / (1.0 - m))
    elseif (m == 0.5) then
      psi = zeta
      psi_int = 0.5 * zeta * zeta
    else
      psi = 1.0 - (1.0 - zeta)**a
      psi_int = zeta - m + m * ((1.0 - zeta)**(1 / m))
    endif
    w_layer = psi
    A_layer = (D_max - D_min)*psi_int
  endif
end subroutine calc_por_layer

function eta_at_uv(e1, e2, interp, h_neglect) result(eatuv)
  real, intent(in) :: e1, e2 ! Interface heights at the adjacent tracer cells [Z ~> m]
  real, intent(in) :: h_neglect ! Negligible thicknesses, often [Z ~> m]
  integer, intent(in) :: interp ! Interpolation method coded by an integer
  real :: eatuv

  select case (interp)
    case (ETA_INTERP_MAX)   ! The shallower interface height
      eatuv = max(e1, e2)
    case (ETA_INTERP_MIN)   ! The deeper interface height
      eatuv = min(e1, e2)
    case (ETA_INTERP_ARITH) ! Arithmetic mean
      eatuv = 0.5 * (e1 + e2)
    case (ETA_INTERP_HARM)  ! Harmonic mean
      eatuv = 2.0 * e1 * e2 / (e1 + e2 + h_neglect)
    case default
      call MOM_error(FATAL, "porous_widths::eta_at_uv: "//&
                     "invalid value for eta interpolation method.")
  end select
end function eta_at_uv

subroutine porous_barriers_init(Time, US, param_file, diag, CS)
  type(porous_barrier_CS), intent(inout) :: CS !< Module control structure
  type(param_file_type),   intent(in)    :: param_file  !< structure indicating parameter file to parse
  type(time_type),         intent(in)    :: Time !< Current model time
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type

  !> This include declares and sets the variable "version".
# include "version_variable.h"
  ! local variables
  character(len=40) :: mdl = "MOM_porous_barriers"  ! This module's name.
  character(len=20) :: interp_method ! String storing eta interpolation method

  CS%initialized = .true.
  CS%diag => diag

  call log_version(param_file, mdl, version, "", log_to_all=.true., layout=.false., &
                   debugging=.false.)
  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false.)
  call get_param(param_file, mdl, "PORBAR_MASKING_DEPTH", CS%mask_depth, &
                 "The depth below which porous barrier is not applied.  "//&
                 "The effective average depths at the velocity cells are used "//&
                 "to test against this criterion.", units="m", default=0.0, &
                 scale=US%m_to_Z)
  CS%mask_depth = -CS%mask_depth
  call get_param(param_file, mdl, "PORBAR_ETA_INTERP", interp_method, &
                 "A string describing the method that decicdes how the "//&
                 "interface heights at the velocity points are calculated. "//&
                 "Valid values are:\n"//&
                 "\t MAX (the default) - maximum of the adjacent cells \n"//&
                 "\t MIN - minimum of the adjacent cells \n"//&
                 "\t ARITHMETIC - arithmetic mean of the adjacent cells \n"//&
                 "\t HARMOINIC - harmonic mean of the adjacent cells \n", &
                 default=ETA_INTERP_MAX_STRING) ! do_not_log=.not.CS%use_por_bar)
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
