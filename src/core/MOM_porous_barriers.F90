!> Module for calculating curve fit for porous topography.
!written by sjd
module MOM_porous_barriers

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs, porous_barrier_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_interface_heights, only : find_eta
use MOM_time_manager,         only : time_type
use MOM_diag_mediator,        only : register_diag_field, diag_ctrl, post_data
use MOM_file_parser,          only : param_file_type, get_param
use MOM_unit_scaling,         only : unit_scale_type
use MOM_debugging,            only : hchksum, uvchksum

implicit none ; private

public porous_widths, porous_barriers_init

#include <MOM_memory.h>

type, public :: porous_barrier_CS; private
  logical :: initialized = .false.  !< True if this control structure has been initialized.
  type(diag_ctrl), pointer :: diag => Null()   !< A structure to regulate diagnostic output timing
  logical :: debug  !< If true, write verbose checksums for debugging purposes.
  real :: mask_depth  !< The depth below which porous barrier is not applied.
  integer :: id_por_layer_widthU = -1, id_por_layer_widthV = -1, id_por_face_areaU = -1, id_por_face_areaV = -1
end type porous_barrier_CS

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
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1):: eta    !< layer interface heights [Z ~> m or 1/eta_to_m].
  integer :: i, j, k, nk, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  real :: dmask
  real :: w_layer, & ! fractional open width of layer interface [nondim]
          A_layer, & ! integral of fractional open width from bottom to current layer[Z ~> m]
          A_layer_prev, & ! integral of fractional open width from bottom to previous layer [Z ~> m]
          eta_s, & ! layer height used for fit [Z ~> m]
          eta_prev ! interface height of previous layer [Z ~> m]

  if (.not.CS%initialized) call MOM_error(FATAL, &
      "MOM_Porous_barrier: Module must be initialized before it is used.")

  is = G%isc; ie = G%iec; js = G%jsc; je = G%jec; nk = GV%ke
  Isq = G%IscB; Ieq = G%IecB; Jsq = G%JscB; Jeq = G%JecB

  dmask = CS%mask_depth

  !eta is zero at surface and decreases downward
  !currently no treatment for using optional find_eta arguments if present
  call find_eta(h, tv, G, GV, US, eta, halo_size=1)

  do j=js,je; do I=Isq,Ieq
    if (G%porous_DavgU(I,j) < dmask) then
      do K = nk+1,1,-1
        eta_s = max(eta(i,j,K), eta(i+1,j,K)) !take shallower layer height
        if (eta_s <= G%porous_DminU(I,j)) then
          pbv%por_layer_widthU(I,j,K) = 0.0
          A_layer_prev = 0.0
          if (K < nk+1) then
            pbv%por_face_areaU(I,j,k) = 0.0; endif
        else
          call calc_por_layer(G%porous_DminU(I,j), G%porous_DmaxU(I,j), &
            G%porous_DavgU(I,j), eta_s, w_layer, A_layer)
          pbv%por_layer_widthU(I,j,K) = w_layer
          if (k <= nk) then
            if ((eta_s - eta_prev) > 0.0) then
              pbv%por_face_areaU(I,j,k) = (A_layer - A_layer_prev)/&
                   (eta_s-eta_prev)
            else
              pbv%por_face_areaU(I,j,k) = 0.0; endif
          endif
          eta_prev = eta_s
          A_layer_prev = A_layer
        endif
      enddo
    endif
  enddo; enddo

  do J=Jsq,Jeq; do i=is,ie
    if (G%porous_DavgV(i,J) < dmask) then
      do K = nk+1,1,-1
        eta_s = max(eta(i,j,K), eta(i,j+1,K)) !take shallower layer height
        if (eta_s <= G%porous_DminV(i,J)) then
          pbv%por_layer_widthV(i,J,K) = 0.0
          A_layer_prev = 0.0
          if (K < nk+1) then
            pbv%por_face_areaV(i,J,k) = 0.0; endif
        else
          call calc_por_layer(G%porous_DminV(i,J), G%porous_DmaxV(i,J), &
            G%porous_DavgV(i,J), eta_s, w_layer, A_layer)
          pbv%por_layer_widthV(i,J,K) = w_layer
          if (k <= nk) then
            if ((eta_s - eta_prev) > 0.0) then
              pbv%por_face_areaV(i,J,k) = (A_layer - A_layer_prev)/&
                   (eta_s-eta_prev)
            else
              pbv%por_face_areaU(I,j,k) = 0.0; endif
          endif
          eta_prev = eta_s
          A_layer_prev = A_layer
        endif
      enddo
    endif
  enddo; enddo

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
end subroutine porous_widths

!> subroutine to calculate the profile fit for a single layer in a column
subroutine calc_por_layer(D_min, D_max, D_avg, eta_layer, w_layer, A_layer)

  real,            intent(in)  :: D_min !< minimum topographic height [Z ~> m]
  real,            intent(in)  :: D_max !< maximum topographic height [Z ~> m]
  real,            intent(in)  :: D_avg !< mean topographic height [Z ~> m]
  real,            intent(in)  :: eta_layer !< height of interface [Z ~> m]
  real,            intent(out) :: w_layer !< frac. open interface width of current layer [nondim]
  real,            intent(out) :: A_layer !< frac. open face area of current layer [Z ~> m]
  !local variables
  real m, a, &             !convenience constant for fit [nondim]
       zeta, &             !normalized vertical coordinate [nondim]
       psi, &              !fractional width of layer between D_min and D_max [nondim]
       psi_int             !integral of psi from 0 to zeta

  !three parameter fit from Adcroft 2013
  m = (D_avg - D_min)/(D_max - D_min)
  a = (1. - m)/m

  zeta = (eta_layer - D_min)/(D_max - D_min)

  if (eta_layer <= D_min) then
    w_layer = 0.0
    A_layer = 0.0
  elseif (eta_layer >= D_max) then
    w_layer = 1.0
    A_layer = eta_layer - D_avg
  else
    if (m < 0.5) then
      psi = zeta**(1./a)
      psi_int = (1.-m)*zeta**(1./(1.-m))
    elseif (m == 0.5) then
      psi = zeta
      psi_int = 0.5*zeta*zeta
    else
      psi = 1. - (1. - zeta)**a
      psi_int = zeta - m + m*((1-zeta)**(1/m))
    endif
    w_layer = psi
    A_layer = (D_max - D_min)*psi_int
  endif


end subroutine calc_por_layer

subroutine porous_barriers_init(Time, US, param_file, diag, CS)
  type(porous_barrier_CS), intent(inout) :: CS
  type(param_file_type),   intent(in)    :: param_file  !< structure indicating parameter file to parse
  type(time_type),         intent(in)    :: Time !< Current model time
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(unit_scale_type),   intent(in)    :: US

  character(len=40)  :: mdl = "MOM_porous_barriers"  ! This module's name.
  CS%initialized = .true.
  CS%diag => diag

  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false.)
  call get_param(param_file, mdl, "POROUS_BARRIER_MASKING_DEPTH", CS%mask_depth,  &
                 "The depth below which porous barrier is not applied.  "//&
                 "This criterion is tested against TOPO_AT_VEL_VARNAME_U_AVE and TOPO_AT_VEL_VARNAME_V_AVE.", &
                 units="m", default=0.0, scale=US%m_to_Z)

  CS%id_por_layer_widthU = register_diag_field('ocean_model', 'por_layer_widthU', diag%axesCui, Time, &
     'Porous barrier open width fraction (at the layer interfaces) of the u-faces', 'nondim')
  CS%id_por_layer_widthV = register_diag_field('ocean_model', 'por_layer_widthV', diag%axesCvi, Time, &
     'Porous barrier open width fraction (at the layer interfaces) of the v-faces', 'nondim')
  CS%id_por_face_areaU = register_diag_field('ocean_model', 'por_face_areaU', diag%axesCuL, Time, &
     'Porous barrier open area fraction (layer averaged) of U-faces', 'nondim')
  CS%id_por_face_areaV = register_diag_field('ocean_model', 'por_face_areaV', diag%axesCvL, Time, &
     'Porous barrier open area fraction (layer averaged) of V-faces', 'nondim')
end subroutine

end module MOM_porous_barriers
