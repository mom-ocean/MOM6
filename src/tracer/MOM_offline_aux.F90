!> Contains routines related to offline transport of tracers. These routines are likely to be called from
!> the MOM_offline_main module
module MOM_offline_aux

! This file is part of MOM6. See LICENSE.md for the license.

use data_override_mod,    only : data_override_init, data_override
use MOM_time_manager,     only : time_type, operator(-)
use MOM_debugging,        only : check_column_integrals
use MOM_domains,          only : pass_var, pass_vector, To_All
use MOM_diag_vkernels,    only : reintegrate_column
use MOM_error_handler,    only : callTree_enter, callTree_leave, MOM_error, FATAL, WARNING, is_root_pe
use MOM_grid,             only : ocean_grid_type
use MOM_io,               only : MOM_read_data, MOM_read_vector, CENTER
use MOM_verticalGrid,     only : verticalGrid_type
use MOM_file_parser,      only : get_param, log_version, param_file_type
use astronomy_mod,        only : orbital_time, diurnal_solar, daily_mean_solar
use MOM_variables,        only : vertvisc_type
use MOM_forcing_type,     only : forcing
use MOM_opacity,          only : optics_type
use MOM_diag_mediator,    only : post_data
use MOM_forcing_type,     only : forcing

implicit none ; private

public update_offline_from_files
public update_offline_from_arrays
public update_h_horizontal_flux
public update_h_vertical_flux
public limit_mass_flux_3d
public distribute_residual_uh_barotropic
public distribute_residual_vh_barotropic
public distribute_residual_uh_upwards
public distribute_residual_vh_upwards
public next_modulo_time
public offline_add_diurnal_sw

#include "MOM_memory.h"
#include "version_variable.h"

contains

!> This updates thickness based on the convergence of horizontal mass fluxes
!! NOTE: Only used in non-ALE mode
subroutine update_h_horizontal_flux(G, GV, uhtr, vhtr, h_pre, h_new)
  type(ocean_grid_type),   pointer       :: G     !< ocean grid structure
  type(verticalGrid_type), pointer       :: GV    !< ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: uhtr  !< Accumulated mass flux through zonal face [kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(in)    :: vhtr  !< Accumulated mass flux through meridional face [kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: h_pre !< Previous layer thicknesses [kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: h_new !< Updated layer thicknesses [kg m-2].

  ! Local variables
  integer :: i, j, k, m, is, ie, js, je, nz
  ! Set index-related variables for fields on T-grid
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

  do k = 1, nz
    do i=is-1,ie+1 ; do j=js-1,je+1

      h_new(i,j,k) = max(0.0, G%US%L_to_m**2*G%areaT(i,j)*h_pre(i,j,k) + &
        ((uhtr(I-1,j,k) - uhtr(I,j,k)) + (vhtr(i,J-1,k) - vhtr(i,J,k))))

      ! In the case that the layer is now dramatically thinner than it was previously,
      ! add a bit of mass to avoid truncation errors.  This will lead to
      ! non-conservation of tracers
      h_new(i,j,k) = h_new(i,j,k) + &
        max(GV%Angstrom_H, 1.0e-13*h_new(i,j,k) - G%US%L_to_m**2*G%areaT(i,j)*h_pre(i,j,k))

      ! Convert back to thickness
      h_new(i,j,k) = h_new(i,j,k) / (G%US%L_to_m**2*G%areaT(i,j))

    enddo ; enddo
  enddo
end subroutine update_h_horizontal_flux

!> Updates layer thicknesses due to vertical mass transports
!! NOTE: Only used in non-ALE configuration
subroutine update_h_vertical_flux(G, GV, ea, eb, h_pre, h_new)
  type(ocean_grid_type),   pointer       :: G     !< ocean grid structure
  type(verticalGrid_type), pointer       :: GV    !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: ea    !< Mass of fluid entrained from the layer
                                                  !! above within this timestep [kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: eb    !< Mass of fluid entrained from the layer
                                                  !! below within this timestep [kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: h_pre !< Layer thicknesses at the end of the previous
                                                  !! step [kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: h_new !< Updated layer thicknesses [kg m-2].

  ! Local variables
  integer :: i, j, k, m, is, ie, js, je, nz
  ! Set index-related variables for fields on T-grid
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

  ! Update h_new with convergence of vertical mass transports
  do j=js-1,je+1
    do i=is-1,ie+1

      ! Top layer
      h_new(i,j,1) = max(0.0, h_pre(i,j,1) + (eb(i,j,1) - ea(i,j,2) + ea(i,j,1) ))
      h_new(i,j,1) = h_new(i,j,1) + &
          max(0.0, 1.0e-13*h_new(i,j,1) - h_pre(i,j,1))

      ! Bottom layer
!        h_new(i,j,nz) = h_pre(i,j,nz) + (ea(i,j,nz) - eb(i,j,nz-1)+eb(i,j,nz))
      h_new(i,j,nz) = max(0.0, h_pre(i,j,nz) + (ea(i,j,nz) - eb(i,j,nz-1)+eb(i,j,nz)))
      h_new(i,j,nz) = h_new(i,j,nz) + &
          max(0.0, 1.0e-13*h_new(i,j,nz) - h_pre(i,j,nz))

    enddo

    ! Interior layers
    do k=2,nz-1 ; do i=is-1,ie+1

      h_new(i,j,k) = max(0.0, h_pre(i,j,k) + ((ea(i,j,k) - eb(i,j,k-1)) + &
          (eb(i,j,k) - ea(i,j,k+1))))
      h_new(i,j,k) = h_new(i,j,k) + &
        max(0.0, 1.0e-13*h_new(i,j,k) - h_pre(i,j,k))

    enddo ; enddo

  enddo

end subroutine update_h_vertical_flux

!> This routine limits the mass fluxes so that the a layer cannot be completely depleted.
!! NOTE: Only used in non-ALE mode
subroutine limit_mass_flux_3d(G, GV, uh, vh, ea, eb, h_pre)
  type(ocean_grid_type),   pointer       :: G     !< ocean grid structure
  type(verticalGrid_type), pointer       :: GV    !< ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: uh    !< Mass flux through zonal face [kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(inout) :: vh    !< Mass flux through meridional face [kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: ea    !< Mass of fluid entrained from the layer
                                                  !! above within this timestep [kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: eb    !< Mass of fluid entrained from the layer
                                                  !! below within this timestep [kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: h_pre !< Layer thicknesses at the end of the previous
                                                  !! step [kg m-2].

  ! Local variables
  integer :: i, j, k, m, is, ie, js, je, nz
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: top_flux, bottom_flux
  real :: pos_flux, hvol, h_neglect, scale_factor, max_off_cfl

  max_off_cfl =0.5

  ! In this subroutine, fluxes out of the box are scaled away if they deplete
  ! the layer, note that we define the positive direction as flux out of the box.
  ! Hence, uh(I-1) is multipled by negative one, but uh(I) is not

  ! Set index-related variables for fields on T-grid
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

  ! Calculate top and bottom fluxes from ea and eb. Note the explicit negative signs
  ! to enforce the positive out convention
  k = 1
  do j=js-1,je+1 ; do i=is-1,ie+1
    top_flux(i,j,k) = -ea(i,j,k)
    bottom_flux(i,j,k) = -(eb(i,j,k)-ea(i,j,k+1))
  enddo ; enddo

  do k=2, nz-1 ; do j=js-1,je+1 ; do i=is-1,ie+1
    top_flux(i,j,k) = -(ea(i,j,k)-eb(i,j,k-1))
    bottom_flux(i,j,k) = -(eb(i,j,k)-ea(i,j,k+1))
  enddo ; enddo ; enddo

  k=nz
  do j=js-1,je+1 ; do i=is-1,ie+1
    top_flux(i,j,k) = -(ea(i,j,k)-eb(i,j,k-1))
    bottom_flux(i,j,k) = -eb(i,j,k)
  enddo ; enddo


  ! Calculate sum of positive fluxes (negatives applied to enforce convention)
  ! in a given cell and scale it back if it would deplete a layer
  do k = 1, nz ; do j=js-1,je+1 ; do i=is-1,ie+1

    hvol = h_pre(i,j,k)*G%US%L_to_m**2*G%areaT(i,j)
    pos_flux  = max(0.0,-uh(I-1,j,k)) + max(0.0, -vh(i,J-1,k)) + &
      max(0.0, uh(I,j,k)) + max(0.0, vh(i,J,k)) + &
      max(0.0, top_flux(i,j,k)*G%US%L_to_m**2*G%areaT(i,j)) + max(0.0, bottom_flux(i,j,k)*G%US%L_to_m**2*G%areaT(i,j))

    if (pos_flux>hvol .and. pos_flux>0.0) then
      scale_factor = ( hvol )/pos_flux*max_off_cfl
    else ! Don't scale
      scale_factor = 1.0
    endif

    ! Scale horizontal fluxes
    if (-uh(I-1,j,k)>0) uh(I-1,j,k) = uh(I-1,j,k)*scale_factor
    if (uh(I,j,k)>0)    uh(I,j,k)   = uh(I,j,k)*scale_factor
    if (-vh(i,J-1,k)>0) vh(i,J-1,k) = vh(i,J-1,k)*scale_factor
    if (vh(i,J,k)>0)    vh(i,J,k)   = vh(i,J,k)*scale_factor

    if (k>1 .and. k<nz) then
    ! Scale interior layers
      if (top_flux(i,j,k)>0.0) then
        ea(i,j,k) = ea(i,j,k)*scale_factor
        eb(i,j,k-1) = eb(i,j,k-1)*scale_factor
      endif
      if (bottom_flux(i,j,k)>0.0) then
        eb(i,j,k) = eb(i,j,k)*scale_factor
        ea(i,j,k+1) = ea(i,j,k+1)*scale_factor
      endif
    ! Scale top layer
    elseif (k==1) then
      if (top_flux(i,j,k)>0.0)    ea(i,j,k) = ea(i,j,k)*scale_factor
      if (bottom_flux(i,j,k)>0.0) then
        eb(i,j,k)   = eb(i,j,k)*scale_factor
        ea(i,j,k+1) = ea(i,j,k+1)*scale_factor
      endif
    ! Scale bottom layer
    elseif (k==nz) then
      if (top_flux(i,j,k)>0.0) then
        ea(i,j,k)   = ea(i,j,k)*scale_factor
        eb(i,j,k-1) = eb(i,j,k-1)*scale_factor
      endif
      if (bottom_flux(i,j,k)>0.0) eb(i,j,k)=eb(i,j,k)*scale_factor
    endif
  enddo ; enddo ; enddo

end subroutine limit_mass_flux_3d

!> In the case where offline advection has failed to converge, redistribute the u-flux
!! into remainder of the water column as a barotropic equivalent
subroutine distribute_residual_uh_barotropic(G, GV, hvol, uh)
  type(ocean_grid_type),   pointer       :: G    !< ocean grid structure
  type(verticalGrid_type), pointer       :: GV   !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                           intent(in   ) :: hvol !< Mass of water in the cells at the end
                                                 !! of the previous timestep [kg]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: uh   !< Zonal mass transport within a timestep [kg]

  real, dimension(SZIB_(G),SZK_(G))   :: uh2d
  real, dimension(SZIB_(G))           :: uh2d_sum
  real, dimension(SZI_(G),SZK_(G))    :: h2d
  real, dimension(SZI_(G))            :: h2d_sum

  integer :: i, j, k, m, is, ie, js, je, nz
  real :: uh_neglect

  ! Set index-related variables for fields on T-grid
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

  do j=js,je
    uh2d_sum(:) = 0.0
    ! Copy over uh to a working array and sum up the remaining fluxes in a column
    do k=1,nz ; do i=is-1,ie
      uh2d(I,k) = uh(I,j,k)
      uh2d_sum(I) = uh2d_sum(I) + uh2d(I,k)
    enddo ; enddo

    ! Copy over h to a working array and calculate total column volume
    h2d_sum(:) = 0.0
    do k=1,nz ; do i=is-1,ie+1
      h2d(i,k) = hvol(i,j,k)
      if (hvol(i,j,k)>0.) then
        h2d_sum(i) = h2d_sum(i) + h2d(i,k)
      else
        h2d(i,k) = GV%H_subroundoff
      endif
    enddo ; enddo

    ! Distribute flux. Note min/max is intended to make sure that the mass transport
    ! does not deplete a cell
    do i=is-1,ie
      if ( uh2d_sum(I)>0.0 ) then
        do k=1,nz
          uh2d(I,k) = uh2d_sum(I)*(h2d(i,k)/h2d_sum(i))
        enddo
      elseif (uh2d_sum(I)<0.0) then
        do k=1,nz
          uh2d(I,k) = uh2d_sum(I)*(h2d(i+1,k)/h2d_sum(i+1))
        enddo
      else
        do k=1,nz
          uh2d(I,k) = 0.0
        enddo
      endif
      ! Calculate and check that column integrated transports match the original to
      ! within the tolerance limit
      uh_neglect = GV%Angstrom_H*G%US%L_to_m**2 * min(G%areaT(i,j), G%areaT(i+1,j))
      if ( abs(sum(uh2d(I,:))-uh2d_sum(I)) > uh_neglect) &
        call MOM_error(WARNING,"Column integral of uh does not match after "//&
        "barotropic redistribution")
    enddo

    do k=1,nz ; do i=is-1,ie
      uh(I,j,k) = uh2d(I,k)
    enddo ; enddo
  enddo

end subroutine distribute_residual_uh_barotropic

!> Redistribute the v-flux as a barotropic equivalent
subroutine distribute_residual_vh_barotropic(G, GV, hvol, vh)
  type(ocean_grid_type),   pointer       :: G    !< ocean grid structure
  type(verticalGrid_type), pointer       :: GV   !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                           intent(in   ) :: hvol !< Mass of water in the cells at the end
                                                 !! of the previous timestep [kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(inout) :: vh   !< Meridional mass transport within a timestep [kg]

  real, dimension(SZJB_(G),SZK_(G))   :: vh2d
  real, dimension(SZJB_(G))           :: vh2d_sum
  real, dimension(SZJ_(G),SZK_(G))    :: h2d
  real, dimension(SZJ_(G))            :: h2d_sum

  integer :: i, j, k, m, is, ie, js, je, nz
  real :: vh_neglect

  ! Set index-related variables for fields on T-grid
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

  do i=is,ie
    vh2d_sum(:) = 0.0
    ! Copy over uh to a working array and sum up the remaining fluxes in a column
    do k=1,nz ; do j=js-1,je
      vh2d(J,k) = vh(i,J,k)
      vh2d_sum(J) = vh2d_sum(J) + vh2d(J,k)
    enddo ; enddo

    ! Copy over h to a working array and calculate column volume
    h2d_sum(:) = 0.0
    do k=1,nz ; do j=js-1,je+1
      h2d(j,k) = hvol(i,j,k)
      if (hvol(i,j,k)>0.) then
        h2d_sum(j) = h2d_sum(j) + h2d(j,k)
      else
        h2d(j,k) = GV%H_subroundoff
      endif
    enddo ; enddo

    ! Distribute flux evenly throughout a column
    do j=js-1,je
      if ( vh2d_sum(J)>0.0 ) then
        do k=1,nz
          vh2d(J,k) = vh2d_sum(J)*(h2d(j,k)/h2d_sum(j))
        enddo
      elseif (vh2d_sum(J)<0.0) then
        do k=1,nz
          vh2d(J,k) = vh2d_sum(J)*(h2d(j+1,k)/h2d_sum(j+1))
        enddo
      else
        do k=1,nz
          vh2d(J,k) = 0.0
        enddo
      endif
      ! Calculate and check that column integrated transports match the original to
      ! within the tolerance limit
      vh_neglect = GV%Angstrom_H*G%US%L_to_m**2 * min(G%areaT(i,j), G%areaT(i,j+1))
      if ( abs(sum(vh2d(J,:))-vh2d_sum(J)) > vh_neglect) then
          call MOM_error(WARNING,"Column integral of vh does not match after "//&
          "barotropic redistribution")
      endif

    enddo

    do k=1,nz ; do j=js-1,je
      vh(i,J,k) = vh2d(J,k)
    enddo ; enddo
  enddo

end subroutine distribute_residual_vh_barotropic

!> In the case where offline advection has failed to converge, redistribute the u-flux
!! into layers above
subroutine distribute_residual_uh_upwards(G, GV, hvol, uh)
  type(ocean_grid_type),   pointer       :: G     !< ocean grid structure
  type(verticalGrid_type), pointer       :: GV    !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                           intent(in   ) :: hvol  !< Mass of water in the cells at the end
                                                  !! of the previous timestep [kg]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: uh    !< Zonal mass transport within a timestep [kg]

  real, dimension(SZIB_(G),SZK_(G))   :: uh2d
  real, dimension(SZI_(G),SZK_(G))    :: h2d

  real  :: uh_neglect, uh_remain, uh_add, uh_sum, uh_col, uh_max
  real  :: hup, hdown, hlos, min_h
  integer :: i, j, k, m, is, ie, js, je, nz, k_rev

  ! Set index-related variables for fields on T-grid
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

  min_h = GV%Angstrom_H*0.1

  do j=js,je
    ! Copy over uh and cell volume to working arrays
    do k=1,nz ; do i=is-2,ie+1
      uh2d(I,k) = uh(I,j,k)
    enddo ; enddo
    do k=1,nz ; do i=is-1,ie+1
      ! Subtract just a little bit of thickness to avoid roundoff errors
      h2d(i,k) = hvol(i,j,k)-min_h*G%US%L_to_m**2*G%areaT(i,j)
    enddo ; enddo

    do i=is-1,ie
      uh_col = SUM(uh2d(I,:)) ! Store original column-integrated transport
      do k=1,nz
        uh_remain = uh2d(I,k)
        uh2d(I,k) = 0.0
        if (abs(uh_remain)>0.0) then
          do k_rev = k,1,-1
            uh_sum = uh_remain + uh2d(I,k_rev)
            if (uh_sum<0.0) then ! Transport to the left
              hup = h2d(i+1,k_rev)
              hlos = max(0.0,uh2d(I+1,k_rev))
              if ((((hup - hlos) + uh_sum) < 0.0) .and. &
                  ((0.5*hup + uh_sum) < 0.0)) then
                uh2d(I,k_rev) = min(-0.5*hup,-hup+hlos,0.0)
                uh_remain = uh_sum - uh2d(I,k_rev)
              else
                uh2d(I,k_rev) = uh_sum
                uh_remain = 0.0
                exit
              endif
            else ! Transport to the right
              hup = h2d(i,k_rev)
              hlos = max(0.0,-uh2d(I-1,k_rev))
              if ((((hup - hlos) - uh_sum) < 0.0) .and. &
                  ((0.5*hup - uh_sum) < 0.0)) then
                uh2d(I,k_rev) = max(0.5*hup,hup-hlos,0.0)
                uh_remain = uh_sum - uh2d(I,k_rev)
              else
                uh2d(I,k_rev) = uh_sum
                uh_remain = 0.0
                exit
              endif
            endif
          enddo ! k_rev
        endif

        if (abs(uh_remain)>0.0) then
          if (k<nz) then
            uh2d(I,k+1) = uh2d(I,k+1) + uh_remain
          else
            uh2d(I,k) = uh2d(I,k) + uh_remain
            call MOM_error(WARNING,"Water column cannot accommodate UH redistribution. Tracer may not be conserved")
          endif
        endif
      enddo ! k-loop

      ! Calculate and check that column integrated transports match the original to
      ! within the tolerance limit
      uh_neglect = GV%Angstrom_H*G%US%L_to_m**2 * min(G%areaT(i,j), G%areaT(i+1,j))
      if (abs(uh_col - sum(uh2d(I,:)))>uh_neglect) then
        call MOM_error(WARNING,"Column integral of uh does not match after "//&
        "upwards redistribution")
      endif

    enddo ! i-loop

    do k=1,nz ; do i=is-1,ie
      uh(I,j,k) = uh2d(I,k)
    enddo ; enddo
  enddo

end subroutine distribute_residual_uh_upwards

!> In the case where offline advection has failed to converge, redistribute the u-flux
!! into layers above
subroutine distribute_residual_vh_upwards(G, GV, hvol, vh)
  type(ocean_grid_type),   pointer       :: G     !< ocean grid structure
  type(verticalGrid_type), pointer       :: GV    !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                           intent(in   ) :: hvol  !< Mass of water in the cells at the end
                                                  !! of the previous timestep [kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(inout) :: vh    !< Meridional mass transport within a timestep [kg]

  real, dimension(SZJB_(G),SZK_(G))   :: vh2d
  real, dimension(SZJB_(G))           :: vh2d_sum
  real, dimension(SZJ_(G),SZK_(G))    :: h2d
  real, dimension(SZJ_(G))            :: h2d_sum

  real  :: vh_neglect, vh_remain, vh_col, vh_sum
  real  :: hup, hlos, min_h
  integer :: i, j, k, m, is, ie, js, je, nz, k_rev

  ! Set index-related variables for fields on T-grid
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

  min_h = 0.1*GV%Angstrom_H

  do i=is,ie
    ! Copy over uh and cell volume to working arrays
    do k=1,nz ; do j=js-2,je+1
      vh2d(J,k) = vh(i,J,k)
    enddo ; enddo
    do k=1,nz ; do j=js-1,je+1
      h2d(j,k) = hvol(i,j,k)-min_h*G%US%L_to_m**2*G%areaT(i,j)
    enddo ; enddo

    do j=js-1,je
      vh_col = SUM(vh2d(J,:))
      do k=1,nz
        vh_remain = vh2d(J,k)
        vh2d(J,k) = 0.0
        if (abs(vh_remain)>0.0) then
          do k_rev = k,1,-1
            vh_sum = vh_remain + vh2d(J,k_rev)
            if (vh_sum<0.0) then ! Transport to the left
              hup = h2d(j+1,k_rev)
              hlos = MAX(0.0,vh2d(J+1,k_rev))
              if ((((hup - hlos) + vh_sum) < 0.0) .and. &
                  ((0.5*hup + vh_sum) < 0.0)) then
                vh2d(J,k_rev) = MIN(-0.5*hup,-hup+hlos,0.0)
                vh_remain = vh_sum - vh2d(J,k_rev)
              else
                vh2d(J,k_rev) = vh_sum
                vh_remain = 0.0
                exit
              endif
            else ! Transport to the right
              hup = h2d(j,k_rev)
              hlos = MAX(0.0,-vh2d(J-1,k_rev))
              if ((((hup - hlos) - vh_sum) < 0.0) .and. &
                  ((0.5*hup - vh_sum) < 0.0)) then
                vh2d(J,k_rev) = MAX(0.5*hup,hup-hlos,0.0)
                vh_remain = vh_sum - vh2d(J,k_rev)
              else
                vh2d(J,k_rev) = vh_sum
                vh_remain = 0.0
                exit
              endif
            endif

          enddo ! k_rev
        endif

        if (abs(vh_remain)>0.0) then
         if (k<nz) then
            vh2d(J,k+1) = vh2d(J,k+1) + vh_remain
          else
            vh2d(J,k) = vh2d(J,k) + vh_remain
            call MOM_error(WARNING,"Water column cannot accommodate VH redistribution. Tracer will not be conserved")
          endif
        endif ! k-loop
      enddo

      ! Calculate and check that column integrated transports match the original to
      ! within the tolerance limit
      vh_neglect = GV%Angstrom_H*G%US%L_to_m**2 * min(G%areaT(i,j), G%areaT(i,j+1))
      if ( ABS(vh_col-SUM(vh2d(J,:))) > vh_neglect) then
        call MOM_error(WARNING,"Column integral of vh does not match after "//&
                               "upwards redistribution")
      endif
    enddo

    do k=1,nz ; do j=js-1,je
      vh(i,J,k) = vh2d(J,k)
    enddo ; enddo
  enddo

end subroutine distribute_residual_vh_upwards

!> add_diurnal_SW adjusts the shortwave fluxes in an forcying_type variable
!! to add a synthetic diurnal cycle. Adapted from SIS2
subroutine offline_add_diurnal_SW(fluxes, G, Time_start, Time_end)
  type(forcing),         intent(inout) :: fluxes !< The type with atmospheric fluxes to be adjusted.
  type(ocean_grid_type), intent(in)    :: G      !< The ocean lateral grid type.
  type(time_type),       intent(in)    :: Time_start !< The start time for this step.
  type(time_type),       intent(in)    :: Time_end   !< The ending time for this step.

  real :: diurnal_factor, time_since_ae, rad
  real :: fracday_dt, fracday_day
  real :: cosz_day, cosz_dt, rrsun_day, rrsun_dt
  type(time_type) :: dt_here

  integer :: i, j, k, i2, j2, isc, iec, jsc, jec, i_off, j_off

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  i_off = LBOUND(fluxes%sens,1) - G%isc ; j_off = LBOUND(fluxes%sens,2) - G%jsc

  !   Orbital_time extracts the time of year relative to the northern
  ! hemisphere autumnal equinox from a time_type variable.
  time_since_ae = orbital_time(Time_start)
  dt_here = Time_end - Time_start
  rad = acos(-1.)/180.

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,rad,Time_start,dt_here,time_since_ae, &
!$OMP                                  fluxes,i_off,j_off) &
!$OMP                          private(i,j,i2,j2,k,cosz_dt,fracday_dt,rrsun_dt, &
!$OMP                                  fracday_day,cosz_day,rrsun_day,diurnal_factor)
  do j=jsc,jec ; do i=isc,iec
!    Per Rick Hemler:
!      Call diurnal_solar with dtime=dt_here to get cosz averaged over dt_here.
!      Call daily_mean_solar to get cosz averaged over a day.  Then
!      diurnal_factor = cosz_dt_ice*fracday_dt_ice*rrsun_dt_ice /
!                       cosz_day*fracday_day*rrsun_day

    call diurnal_solar(G%geoLatT(i,j)*rad, G%geoLonT(i,j)*rad, Time_start, cosz=cosz_dt, &
                       fracday=fracday_dt, rrsun=rrsun_dt, dt_time=dt_here)
    call daily_mean_solar(G%geoLatT(i,j)*rad, time_since_ae, cosz_day, fracday_day, rrsun_day)
    diurnal_factor = cosz_dt*fracday_dt*rrsun_dt / &
                     max(1e-30, cosz_day*fracday_day*rrsun_day)

    i2 = i+i_off ; j2 = j+j_off
    fluxes%sw(i2,j2) = fluxes%sw(i2,j2) * diurnal_factor
    fluxes%sw_vis_dir(i2,j2) = fluxes%sw_vis_dir(i2,j2) * diurnal_factor
    fluxes%sw_vis_dif(i2,j2) = fluxes%sw_vis_dif(i2,j2) * diurnal_factor
    fluxes%sw_nir_dir(i2,j2) = fluxes%sw_nir_dir(i2,j2) * diurnal_factor
    fluxes%sw_nir_dif(i2,j2) = fluxes%sw_nir_dif(i2,j2) * diurnal_factor
  enddo ; enddo

end subroutine offline_add_diurnal_sw

!> Controls the reading in 3d mass fluxes, diffusive fluxes, and other fields stored
!! in a previous integration of the online model
subroutine update_offline_from_files(G, GV, nk_input, mean_file, sum_file, snap_file, surf_file, h_end, &
    uhtr, vhtr, temp_mean, salt_mean, mld, Kd, fluxes, ridx_sum, ridx_snap, read_mld, read_sw, &
    read_ts_uvh, do_ale_in)

  type(ocean_grid_type),   intent(inout) :: G         !< Horizontal grid type
  type(verticalGrid_type), intent(in   ) :: GV        !< Vertical grid type
  integer,                 intent(in   ) :: nk_input  !< Number of levels in input file
  character(len=*),        intent(in   ) :: mean_file !< Name of file with averages fields
  character(len=*),        intent(in   ) :: sum_file  !< Name of file with summed fields
  character(len=*),        intent(in   ) :: snap_file !< Name of file with snapshot fields
  character(len=*),        intent(in   ) :: surf_file !< Name of file with surface fields
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: uhtr      !< Zonal mass fluxes [kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(inout) :: vhtr      !< Meridional mass fluxes [kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                           intent(inout) :: h_end     !< End of timestep layer thickness
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                           intent(inout) :: temp_mean !< Averaged temperature
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                           intent(inout) :: salt_mean !< Averaged salinity
  real, dimension(SZI_(G),SZJ_(G)),          &
                           intent(inout) :: mld       !< Averaged mixed layer depth
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), &
                           intent(inout) :: Kd        !< Diapycnal diffusivities at interfaces
  type(forcing),           intent(inout) :: fluxes    !< Fields with surface fluxes
  integer,                 intent(in   ) :: ridx_sum  !< Read index for sum, mean, and surf files
  integer,                 intent(in   ) :: ridx_snap !< Read index for snapshot file
  logical,                 intent(in   ) :: read_mld  !< True if reading in MLD
  logical,                 intent(in   ) :: read_sw   !< True if reading in radiative fluxes
  logical,                 intent(in   ) :: read_ts_uvh !< True if reading in uh, vh, and h
  logical,       optional, intent(in   ) :: do_ale_in !< True if using ALE algorithms

  logical :: do_ale
  integer :: i, j, k, is, ie, js, je, nz
  real    :: Initer_vert

  do_ale = .false.
  if (present(do_ale_in) ) do_ale = do_ale_in

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  ! Check if reading in UH, VH, and h_end
  if (read_ts_uvh) then
    h_end(:,:,:) = 0.0
    temp_mean(:,:,:) = 0.0
    salt_mean(:,:,:) = 0.0
    uhtr(:,:,:) = 0.0
    vhtr(:,:,:) = 0.0
    ! Time-summed fields
    call MOM_read_vector(sum_file, 'uhtr_sum', 'vhtr_sum', uhtr(:,:,1:nk_input), &
                         vhtr(:,:,1:nk_input), G%Domain, timelevel=ridx_sum)
    call MOM_read_data(snap_file, 'h_end', h_end(:,:,1:nk_input), G%Domain, &
                       timelevel=ridx_snap,position=CENTER)
    call MOM_read_data(mean_file, 'temp', temp_mean(:,:,1:nk_input), G%Domain, &
                       timelevel=ridx_sum,position=CENTER)
    call MOM_read_data(mean_file, 'salt', salt_mean(:,:,1:nk_input), G%Domain, &
                       timelevel=ridx_sum,position=CENTER)
  endif

  do j=js,je ; do i=is,ie
    if (G%mask2dT(i,j)>0.) then
      temp_mean(:,:,nk_input:nz) = temp_mean(i,j,nk_input)
      salt_mean(:,:,nk_input:nz) = salt_mean(i,j,nk_input)
    endif
  enddo ; enddo

  ! Check if reading vertical diffusivities or entrainment fluxes
  call MOM_read_data( mean_file, 'Kd_interface', Kd(:,:,1:nk_input+1), G%Domain, &
                  timelevel=ridx_sum,position=CENTER)

  ! This block makes sure that the fluxes control structure, which may not be used in the solo_driver,
  ! contains netMassIn and netMassOut which is necessary for the applyTracerBoundaryFluxesInOut routine
  if (do_ale) then
    if (.not. associated(fluxes%netMassOut)) then
      allocate(fluxes%netMassOut(G%isd:G%ied,G%jsd:G%jed))
      fluxes%netMassOut(:,:) = 0.0
    endif
    if (.not. associated(fluxes%netMassIn)) then
      allocate(fluxes%netMassIn(G%isd:G%ied,G%jsd:G%jed))
      fluxes%netMassIn(:,:) = 0.0
    endif

    fluxes%netMassOut(:,:) = 0.0
    fluxes%netMassIn(:,:) = 0.0
    call MOM_read_data(surf_file,'massout_flux_sum',fluxes%netMassOut, G%Domain, &
        timelevel=ridx_sum)
    call MOM_read_data(surf_file,'massin_flux_sum', fluxes%netMassIn,  G%Domain, &
        timelevel=ridx_sum)

    do j=js,je ; do i=is,ie
      if (G%mask2dT(i,j)<1.0) then
        fluxes%netMassOut(i,j) = 0.0
        fluxes%netMassIn(i,j) = 0.0
      endif
    enddo ; enddo

  endif

  if (read_mld) then
    call MOM_read_data(surf_file, 'ePBL_h_ML', mld, G%Domain, timelevel=ridx_sum)
  endif

  if (read_sw) then
    ! Shortwave radiation is only needed for offline mode with biogeochemistry but without the coupler.
    ! Need to double check, but set_opacity seems to only need the sum of the diffuse and
    ! direct fluxes in the visible and near-infrared bands. For convenience, we store the
    ! sum of the direct and diffuse fluxes in the 'dir' field and set the 'dif' fields to zero
    call MOM_read_data(mean_file,'sw_vis', fluxes%sw_vis_dir, G%Domain, &
                       timelevel=ridx_sum, scale=G%US%W_m2_to_QRZ_T)
    call MOM_read_data(mean_file,'sw_nir', fluxes%sw_nir_dir, G%Domain, &
                       timelevel=ridx_sum, scale=G%US%W_m2_to_QRZ_T)
    fluxes%sw_vis_dir(:,:) = fluxes%sw_vis_dir(:,:)*0.5
    fluxes%sw_vis_dif(:,:) = fluxes%sw_vis_dir(:,:)
    fluxes%sw_nir_dir(:,:) = fluxes%sw_nir_dir(:,:)*0.5
    fluxes%sw_nir_dif(:,:) = fluxes%sw_nir_dir(:,:)
    fluxes%sw = (fluxes%sw_vis_dir + fluxes%sw_vis_dif) + (fluxes%sw_nir_dir + fluxes%sw_nir_dif)
    do j=js,je ; do i=is,ie
      if (G%mask2dT(i,j)<1.0) then
        fluxes%sw(i,j) = 0.0
        fluxes%sw_vis_dir(i,j) = 0.0
        fluxes%sw_nir_dir(i,j) = 0.0
        fluxes%sw_vis_dif(i,j) = 0.0
        fluxes%sw_nir_dif(i,j) = 0.0
      endif
    enddo ; enddo
    call pass_var(fluxes%sw,G%Domain)
    call pass_var(fluxes%sw_vis_dir,G%Domain)
    call pass_var(fluxes%sw_vis_dif,G%Domain)
    call pass_var(fluxes%sw_nir_dir,G%Domain)
    call pass_var(fluxes%sw_nir_dif,G%Domain)
  endif

end subroutine update_offline_from_files

!> Fields for offline transport are copied from the stored arrays read during initialization
subroutine update_offline_from_arrays(G, GV, nk_input, ridx_sum, mean_file, sum_file, snap_file, uhtr, vhtr, &
                                      hend, uhtr_all, vhtr_all, hend_all, temp, salt, temp_all, salt_all )
  type(ocean_grid_type),                     intent(inout) :: G         !< Horizontal grid type
  type(verticalGrid_type),                   intent(in   ) :: GV        !< Vertical grid type
  integer,                                   intent(in   ) :: nk_input  !< Number of levels in input file
  integer,                                   intent(in   ) :: ridx_sum  !< Index to read from
  character(len=200),                        intent(in   ) :: mean_file !< Name of file with averages fields
  character(len=200),                        intent(in   ) :: sum_file  !< Name of file with summed fields
  character(len=200),                        intent(in   ) :: snap_file !< Name of file with snapshot fields
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uhtr      !< Zonal mass fluxes [kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vhtr      !< Meridional mass fluxes [kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: hend      !< End of timestep layer thickness [kg m-2]
  real, dimension(:,:,:,:), allocatable,     intent(inout) :: uhtr_all  !< Zonal mass fluxes [kg]
  real, dimension(:,:,:,:), allocatable,     intent(inout) :: vhtr_all  !< Meridional mass fluxes [kg]
  real, dimension(:,:,:,:), allocatable,     intent(inout) :: hend_all  !< End of timestep layer thickness [kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: temp      !< Temperature array
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: salt      !< Salinity array
  real, dimension(:,:,:,:), allocatable,     intent(inout) :: temp_all  !< Temperature array
  real, dimension(:,:,:,:), allocatable,     intent(inout) :: salt_all  !< Salinity array

  integer :: i, j, k, is, ie, js, je, nz
  real, parameter :: fill_value = 0.
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  ! Check that all fields are allocated (this is a redundant check)
  if (.not. allocated(uhtr_all)) &
      call MOM_error(FATAL, "uhtr_all not allocated before call to update_transport_from_arrays")
  if (.not. allocated(vhtr_all)) &
      call MOM_error(FATAL, "vhtr_all not allocated before call to update_transport_from_arrays")
  if (.not. allocated(hend_all)) &
      call MOM_error(FATAL, "hend_all not allocated before call to update_transport_from_arrays")
  if (.not. allocated(temp_all)) &
      call MOM_error(FATAL, "temp_all not allocated before call to update_transport_from_arrays")
  if (.not. allocated(salt_all)) &
      call MOM_error(FATAL, "salt_all not allocated before call to update_transport_from_arrays")

  ! Copy uh, vh, h_end, temp, and salt
  do k=1,nk_input ; do j=js,je ; do i=is,ie
    uhtr(I,j,k) = uhtr_all(I,j,k,ridx_sum)
    vhtr(i,J,k) = vhtr_all(i,J,k,ridx_sum)
    hend(i,j,k) = hend_all(i,j,k,ridx_sum)
    temp(i,j,k) = temp_all(i,j,k,ridx_sum)
    salt(i,j,k) = salt_all(i,j,k,ridx_sum)
  enddo ; enddo ; enddo

  ! Fill the rest of the arrays with 0s (fill_value could probably be changed to a runtime parameter)
  do k=nk_input+1,nz ; do j=js,je ; do i=is,ie
    uhtr(I,j,k) = fill_value
    vhtr(i,J,k) = fill_value
    hend(i,j,k) = fill_value
    temp(i,j,k) = fill_value
    salt(i,j,k) = fill_value
  enddo ; enddo ; enddo

end subroutine update_offline_from_arrays

!> Calculates the next timelevel to read from the input fields. This allows the 'looping'
!! of the fields
function next_modulo_time(inidx, numtime)
  ! Returns the next time interval to be read
  integer                 :: numtime              ! Number of time levels in input fields
  integer                 :: inidx                ! The current time index

  integer                 :: read_index           ! The index in the input files that corresponds
                                                  ! to the current timestep

  integer                 :: next_modulo_time

  read_index = mod(inidx+1,numtime)
  if (read_index < 0)  read_index = inidx-read_index
  if (read_index == 0) read_index = numtime

  next_modulo_time = read_index

end function next_modulo_time

end module MOM_offline_aux

