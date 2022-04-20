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

implicit none ; private

#include <MOM_memory.h>

public porous_widths

!> Calculates curve fit from D_min, D_max, D_avg
interface porous_widths
  module procedure por_widths, calc_por_layer
end interface porous_widths

contains

!> subroutine to assign cell face areas and layer widths for porous topography
subroutine por_widths(h, tv, G, GV, US, eta, pbv, eta_bt, halo_size, eta_to_m)
  !eta_bt, halo_size, eta_to_m not currently used
  !variables needed to call find_eta
  type(ocean_grid_type),                      intent(in)  :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)  :: tv     !< A structure pointing to various
                                                                    !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(out) :: eta    !< layer interface heights
                                                                    !! [Z ~> m] or 1/eta_to_m m).
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)  :: eta_bt !< optional barotropic
             !! variable that gives the "correct" free surface height (Boussinesq) or total water
             !! column mass per unit area (non-Boussinesq).  This is used to dilate the layer.
             !! thicknesses when calculating interfaceheights [H ~> m or kg m-2].
  integer,                          optional, intent(in)  :: halo_size !< width of halo points on
                                                                       !! which to calculate eta.

  real,                             optional, intent(in)  :: eta_to_m  !< The conversion factor from
             !! the units of eta to m; by default this is US%Z_to_m.
  type(porous_barrier_ptrs),           intent(inout) :: pbv  !< porous barrier fractional cell metrics

  !local variables
  integer i, j, k, nk, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real w_layer, & ! fractional open width of layer interface [nondim]
       A_layer, & ! integral of fractional open width from bottom to current layer[Z ~> m]
       A_layer_prev, & ! integral of fractional open width from bottom to previous layer [Z ~> m]
       eta_s, & ! layer height used for fit [Z ~> m]
       eta_prev ! interface height of previous layer [Z ~> m]
  isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed
  IsdB = G%IsdB; IedB = G%IedB; JsdB = G%JsdB; JedB = G%JedB

  !eta is zero at surface and decreases downward

  nk = SZK_(G)

  !currently no treatment for using optional find_eta arguments if present
  call find_eta(h, tv, G, GV, US, eta)

  do j=jsd,jed; do I=IsdB,IedB
    if (G%porous_DavgU(I,j) < 0.) then
      do K = nk+1,1,-1
        eta_s = max(eta(I,j,K), eta(I+1,j,K)) !take shallower layer height
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

  do J=JsdB,JedB; do i=isd,ied
    if (G%porous_DavgV(i,J) < 0.) then
      do K = nk+1,1,-1
        eta_s = max(eta(i,J,K), eta(i,J+1,K)) !take shallower layer height
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

end subroutine por_widths

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

end module MOM_porous_barriers
