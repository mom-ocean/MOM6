!> Contains routines related to offline transport of tracers
module MOM_offline_aux
! This file is part of MOM6. See LICENSE.md for the license.

use data_override_mod,    only : data_override_init, data_override
use MOM_time_manager,     only : time_type, operator(-)
use MOM_domains,          only : pass_var, pass_vector, To_All
use MOM_error_handler,    only : callTree_enter, callTree_leave, MOM_error, FATAL, WARNING, is_root_pe
use MOM_grid,             only : ocean_grid_type
use MOM_verticalGrid,     only : verticalGrid_type
use MOM_file_parser,      only : get_param, log_version, param_file_type
use astronomy_mod,        only : orbital_time, diurnal_solar, daily_mean_solar
use MOM_variables,        only : vertvisc_type
use MOM_forcing_type,     only : forcing
use MOM_shortwave_abs,    only : optics_type
use MOM_diag_mediator,    only : post_data
use MOM_forcing_type,     only : forcing

implicit none

public update_h_horizontal_flux
public update_h_vertical_flux
public limit_mass_flux_3d
public distribute_residual_uh_barotropic
public distribute_residual_vh_barotropic
public distribute_residual_uh_upwards
public distribute_residual_vh_upwards
public offline_add_diurnal_sw

#include "MOM_memory.h"
#include "version_variable.h"

contains

!> This updates thickness based on the convergence of horizontal mass fluxes
!! NOTE: Only used in non-ALE mode
subroutine update_h_horizontal_flux(G, GV, uhtr, vhtr, h_pre, h_new)
  type(ocean_grid_type),    pointer                           :: G
  type(verticalGrid_type),  pointer                           :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)       :: uhtr
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)       :: vhtr
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(in)       :: h_pre
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: h_new

  ! Local variables
  integer :: i, j, k, m, is, ie, js, je, nz
  ! Set index-related variables for fields on T-grid
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

  do k = 1, nz
    do i=is-1,ie+1 ; do j=js-1,je+1

      h_new(i,j,k) = max(0.0, G%areaT(i,j)*h_pre(i,j,k) + &
        ((uhtr(I-1,j,k) - uhtr(I,j,k)) + (vhtr(i,J-1,k) - vhtr(i,J,k))))

      ! In the case that the layer is now dramatically thinner than it was previously,
      ! add a bit of mass to avoid truncation errors.  This will lead to
      ! non-conservation of tracers
      h_new(i,j,k) = h_new(i,j,k) + &
        max(GV%Angstrom, 1.0e-13*h_new(i,j,k) - G%areaT(i,j)*h_pre(i,j,k))

      ! Convert back to thickness
      h_new(i,j,k) = h_new(i,j,k)/G%areaT(i,j)

    enddo ; enddo
  enddo
end subroutine update_h_horizontal_flux

!> Updates layer thicknesses due to vertical mass transports
!! NOTE: Only used in non-ALE configuration
subroutine update_h_vertical_flux(G, GV, ea, eb, h_pre, h_new)
  type(ocean_grid_type),    pointer                           :: G
  type(verticalGrid_type),  pointer                           :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(in)       :: ea
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(in)       :: eb
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(in)       :: h_pre
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: h_new

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
subroutine limit_mass_flux_3d(G, GV, uh, vh, ea, eb, h_pre, max_off_cfl)
  type(ocean_grid_type),    pointer                           :: G
  type(verticalGrid_type),  pointer                           :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: uh
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)    :: vh
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: ea
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: eb
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(in)       :: h_pre
  real,                                      intent(in)       :: max_off_cfl

  ! Local variables
  integer :: i, j, k, m, is, ie, js, je, nz
  real, dimension(SZI_(G),SZJ_(G),SZK_(G))                    :: top_flux, bottom_flux
  real                                                        :: pos_flux, hvol, h_neglect, scale_factor


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

    hvol = h_pre(i,j,k)*G%areaT(i,j)
    pos_flux  = max(0.0,-uh(I-1,j,k)) + max(0.0, -vh(i,J-1,k)) + &
      max(0.0, uh(I,j,k)) + max(0.0, vh(i,J,k)) + &
      max(0.0, top_flux(i,j,k)*G%areaT(i,j)) + max(0.0, bottom_flux(i,j,k)*G%areaT(i,j))

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
      if(top_flux(i,j,k)>0.0) then
        ea(i,j,k) = ea(i,j,k)*scale_factor
        eb(i,j,k-1) = eb(i,j,k-1)*scale_factor
      endif
      if(bottom_flux(i,j,k)>0.0) then
        eb(i,j,k) = eb(i,j,k)*scale_factor
        ea(i,j,k+1) = ea(i,j,k+1)*scale_factor
      endif
    ! Scale top layer
    elseif (k==1) then
      if(top_flux(i,j,k)>0.0)    ea(i,j,k) = ea(i,j,k)*scale_factor
      if(bottom_flux(i,j,k)>0.0) then
        eb(i,j,k)   = eb(i,j,k)*scale_factor
        ea(i,j,k+1) = ea(i,j,k+1)*scale_factor
      endif
    ! Scale bottom layer
    elseif (k==nz) then
      if(top_flux(i,j,k)>0.0) then
        ea(i,j,k)   = ea(i,j,k)*scale_factor
        eb(i,j,k-1) = eb(i,j,k-1)*scale_factor
      endif
      if (bottom_flux(i,j,k)>0.0) eb(i,j,k)=eb(i,j,k)*scale_factor
    endif
  enddo ; enddo ; enddo

end subroutine limit_mass_flux_3d

!> In the case where offline advection has failed to converge, redistribute the u-flux
!! into remainder of the water column as a barotropic equivalent
subroutine distribute_residual_uh_barotropic(G, GV, h, uh)
  type(ocean_grid_type),    pointer                           :: G
  type(verticalGrid_type),  pointer                           :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout)    :: h
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: uh

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

    ! Copy over h to a working array and calculate column height
    h2d_sum(:) = 0.0
    do k=1,nz ; do i=is-2,ie+1
      h2d(i,k) = h(i,j,k)*G%areaT(i,j)
      if(h(i,j,k)>GV%Angstrom) then
        h2d_sum(i) = h2d_sum(i) + h2d(i,k)
      else
        h2d(i,k) = 0.0
      endif
    enddo; enddo;


    ! Distribute flux. Note min/max is intended to make sure that the mass transport
    ! does not deplete a cell
    do i=is-1,ie
      if( uh2d_sum(I)>0.0 ) then
        do k=1,nz
          uh2d(I,k) = min(uh2d_sum(I)*(h2d(i,k)/h2d_sum(i)),h2d(i,k))
        enddo
      elseif (uh2d_sum(I)<0.0) then
        do k=1,nz
          uh2d(I,k) = max(uh2d_sum(I)*(h2d(i+1,k)/h2d_sum(i+1)),-h2d(i+1,k))
        enddo
      else
        do k=1,nz
          uh2d(I,k) = 0.0
        enddo
      endif
      ! Calculate and check that column integrated transports match the original to
      ! within the tolerance limit
      uh_neglect = nz*GV%Angstrom*min(G%areaT(i,j),G%areaT(i+1,j))
      if( abs(sum(uh2d(I,:))-uh2d_sum(I)) > uh_neglect) &
        call MOM_error(WARNING,"Column integral of uh does not match after "//&
        "barotropic redistribution")
    enddo

    ! Update layer thicknesses at the end
    do k=1,nz ; do i=is,ie
      h(i,j,k) = h(i,j,k) + (uh2d(I-1,k) - uh2d(I,k))/G%areaT(i,j)
    enddo ; enddo
    do k=1,nz ; do i=is-1,ie
      uh(I,j,k) = uh2d(I,k)
    enddo ; enddo
  enddo

end subroutine distribute_residual_uh_barotropic

!> Redistribute the v-flux as a barotropic equivalent
subroutine distribute_residual_vh_barotropic(G, GV, h, vh)
  type(ocean_grid_type),    pointer                           :: G
  type(verticalGrid_type),  pointer                           :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout)    :: h
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)    :: vh

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
    do k=1,nz ; do j=js-2,je+1
      h2d(j,k) = h(i,j,k)*G%areaT(i,j)
      if(h(i,j,k)>GV%Angstrom) then
        h2d_sum(j) = h2d_sum(j) + h2d(j,k)
      else
        h2d(j,k) = 0.0
      endif
    enddo; enddo;


    ! Distribute flux. Note min/max is intended to make sure that the mass transport
    ! does not deplete a cell. If this limit is hit for some reason, tracer will
    ! not be conserved
    do j=js-1,je
      if( vh2d_sum(J)>0.0 ) then
        do k=1,nz
          vh2d(J,k) = min(vh2d_sum(J)*(h2d(j,k)/h2d_sum(j)),0.5*h2d(j,k))
        enddo
      elseif (vh2d_sum(J)<0.0) then
        do k=1,nz
          vh2d(J,k) = max(vh2d_sum(J)*(h2d(j+1,k)/h2d_sum(j+1)),-0.5*h2d(j+1,k))
        enddo
      else
        do k=1,nz
          vh2d(J,k) = 0.0
        enddo
      endif
      ! Calculate and check that column integrated transports match the original to
      ! within the tolerance limit
      vh_neglect = nz*GV%Angstrom*min(G%areaT(i,j),G%areaT(i,j+1))
      if( abs(sum(vh2d(J,:))-vh2d_sum(J)) > vh_neglect) &
          call MOM_error(WARNING,"Column integral of vh does not match after "//&
          "barotropic redistribution")

    enddo

    ! Update layer thicknesses at the end.
    ! This may not be needed since the limits on the flux are half of the original thickness
    do k=1,nz ; do j=js,je
      h(i,j,k) = h(i,j,k) + (vh2d(J-1,k) - vh2d(J,k))/G%areaT(i,j)
    enddo ; enddo
    do k=1,nz ; do j=js-1,je
      vh(i,J,k) = vh2d(J,k)
    enddo ; enddo
  enddo

end subroutine distribute_residual_vh_barotropic

!> In the case where offline advection has failed to converge, redistribute the u-flux
!! into layers above
subroutine distribute_residual_uh_upwards(G, GV, h, uh)
  type(ocean_grid_type),    pointer                           :: G
  type(verticalGrid_type),  pointer                           :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout)    :: h
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: uh

  real, dimension(SZIB_(G),SZK_(G))   :: uh2d
  real, dimension(SZI_(G),SZK_(G))    :: h2d

  real  :: uh_neglect, uh_remain, uh_LB, uh_UB, uh_add
  real  :: hup, hdown, hlos, min_h
  integer :: i, j, k, m, is, ie, js, je, nz, k_rev

  ! Set index-related variables for fields on T-grid
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

  min_h = 0.1*GV%Angstrom

  do j=js-1,je
    ! Copy over uh and cell volume to working arrays
    do k=1,nz ; do i=is-1,ie
      uh2d(I,k) = uh(I,j,k)
    enddo ; enddo
    do k=1,nz ; do i=is-2,ie+1
      ! Subtract just a little bit of thickness to avoid roundoff errors
      h2d(i,k) = max(h(i,j,k)*G%areaT(i,j)-min_h*G%areaT(i,j),min_h*G%areaT(i,j))
    enddo ; enddo

    do i=is-1,ie
      do k=1,nz
        uh_remain = uh2d(I,k)
        uh_neglect = GV%H_subroundoff*min(G%areaT(i,j),G%areaT(i+1,j))
        if(uh_remain<-uh_neglect) then
          ! Set the mass flux to zero. This will be refilled in the first iteration
          uh2d(I,k) = 0.0
          do k_rev=k,1,-1
            ! This lower bound only allows half of the layer to be depleted
            uh_LB = -0.5*h2d(i+1,k_rev)
            ! You can either add the difference between the lower bound and the
            ! current uh, or the remaining mass transport to be distributed.
            ! The max is there because it represents the minimum of these two with respect
            ! to magnitude. The minimum is to guard against the case where uh2d>uh_LB
            ! not quite the same potentially because of roundoff error
            uh_add = min(max(uh_LB-uh2d(I,k_rev), uh_remain),0.0)
            uh_remain = uh_remain - uh_add
            uh2d(I,k_rev) = uh2d(I,k_rev) + uh_add
            if(uh_remain>-uh_neglect) exit
          enddo
        elseif (uh_remain>uh_neglect) then
          ! Set the amount in the layer with remaining fluxes to zero. This will be reset
          ! in the first iteration of the redistribution loop
          uh2d(I,k) = 0.0
          ! Loop to distribute remaining flux in layers above
          do k_rev=k,1,-1
            ! This lower bound only allows half of the layer to be depleted
            uh_UB = 0.5*h2d(i,k_rev)
            uh_add = max(min(uh_UB-uh2d(I,k_rev), uh_remain), 0.0)
            uh_remain = uh_remain - uh_add
            uh2d(I,k_rev) = uh2d(I,k_rev) + uh_add
            if(uh_remain<uh_neglect) exit
          enddo
          ! Check to see if there's any mass flux left. If so, put it in the layer beneath,
          ! unless we've bottomed out
        endif
        if(abs(uh_remain)>uh_neglect) then
          if(k<nz) then
            uh2d(I,k+1) = uh2d(I,k+1) + uh_remain
          else
            call MOM_error(WARNING,"Water column cannot accommodate UH redistribution. Tracer will not be conserved")
          endif
        endif

      enddo

    enddo

    ! Update layer thicknesses at the end
    do k=1,nz ; do i=is,ie
      h(i,j,k) = (h(i,j,k)*G%areaT(i,j) + (uh2d(I-1,k) - uh2d(I,k)))/G%areaT(i,j)
    enddo ; enddo
    do k=1,nz ; do i=is-1,ie
      uh(I,j,k) = uh2d(I,k)
    enddo ; enddo
  enddo

end subroutine distribute_residual_uh_upwards

!> In the case where offline advection has failed to converge, redistribute the u-flux
!! into layers above
subroutine distribute_residual_vh_upwards(G, GV, h, vh)
  type(ocean_grid_type),    pointer                          :: G
  type(verticalGrid_type),  pointer                          :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout)   :: h
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)   :: vh

  real, dimension(SZJB_(G),SZK_(G))   :: vh2d
  real, dimension(SZJB_(G))           :: vh2d_sum
  real, dimension(SZJ_(G),SZK_(G))    :: h2d
  real, dimension(SZJ_(G))            :: h2d_sum

  real  :: vh_neglect, vh_remain, vh_max, vh_add, vh_UB, vh_LB, vh_sum
  real  :: hup, hlos, min_h
  integer :: i, j, k, m, is, ie, js, je, nz, k_rev

  ! Set index-related variables for fields on T-grid
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

  min_h = 0.1*GV%Angstrom

  do i=is-1,ie
    ! Copy over uh and cell volume to working arrays
    do k=1,nz ; do j=js-1,je
      vh2d(J,k) = vh(i,J,k)
    enddo ; enddo
    do k=1,nz ; do j=js-2,je+1
      h2d(j,k) = (h(i,j,k)-min_h)*G%areaT(i,j)
    enddo ; enddo

    do j=js-1,je
      do k=1,nz
        vh_remain = vh2d(J,k)
        vh_neglect = GV%H_subroundoff*min(G%areaT(i,j),G%areaT(i,j+1))
        if(vh_remain<-vh_neglect) then
          ! Set the mass flux to zero. This will be refilled in the first iteration
          vh2d(J,k) = 0.0
          do k_rev=k,1,-1
            ! This lower bound only allows half of the layer to be depleted
            vh_LB = -0.5*h2d(j+1,k_rev)
            ! You can either add the difference between the lower bound and the
            ! current uh, or the remaining mass transport to be distributed.
            ! The max is there because it represents the minimum of these two with respect
            ! to magnitude. The minimum is to guard against the case where uh2d>uh_LB
            ! not quite the same potentially because of roundoff error
            vh_add = min(max(vh_LB-vh2d(J,k_rev), vh_remain),0.0)
            vh_remain = vh_remain - vh_add
            vh2d(J,k_rev) = vh2d(J,k_rev) + vh_add
            if(vh_remain>-vh_neglect) exit
          enddo
        elseif (vh_remain>vh_neglect) then
          ! Set the amount in the layer with remaining fluxes to zero. This will be reset
          ! in the first iteration of the redistribution loop
          vh2d(J,k) = 0.0
          ! Loop to distribute remaining flux in layers above
          do k_rev=k,1,-1
            ! This lower bound only allows half of the layer to be depleted
            vh_UB = 0.5*h2d(j,k_rev)
            vh_add = max(min(vh_UB-vh2d(J,k_rev), vh_remain), 0.0)
            vh_remain = vh_remain - vh_add
            vh2d(J,k_rev) = vh2d(J,k_rev) + vh_add
            if(vh_remain<vh_neglect) exit
          enddo
          ! Check to see if there's any mass flux left. If so, put it in the layer beneath,
          ! unless we've bottomed out
        endif
        if(abs(vh_remain)>vh_neglect) then
          if(k<nz) then
            vh2d(J,k+1) = vh2d(J,k+1) + vh_remain
          else
            call MOM_error(WARNING,"Water column cannot accommodate UH redistribution. Tracer will not be conserved")
          endif
        endif
      enddo

    enddo

    ! Update layer thicknesses at the end
    do k=1,nz ; do j=js,je
      h(i,j,k) = (h(i,j,k)*G%areaT(i,j) + (vh2d(J-1,k) - vh2d(J,k)))/G%areaT(i,j)
    enddo ; enddo
    do k=1,nz ; do j=js-1,je
      vh(i,J,k) = vh2d(J,k)
    enddo ; enddo
  enddo


end subroutine distribute_residual_vh_upwards

!> add_diurnal_SW adjusts the shortwave fluxes in an forcying_type variable
!! to add a synthetic diurnal cycle. Adapted from SIS2
subroutine offline_add_diurnal_SW(fluxes, G, Time_start, Time_end)
  type(forcing),                 intent(inout) :: fluxes !< The type with atmospheric fluxes to be adjusted.
  type(ocean_grid_type),         intent(in)    :: G   !< The sea-ice lateral grid type.
  type(time_type),               intent(in)    :: Time_start !< The start time for this step.
  type(time_type),               intent(in)    :: Time_end   !< The ending time for this step.

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
    call daily_mean_solar (G%geoLatT(i,j)*rad, time_since_ae, cosz_day, fracday_day, rrsun_day)
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

end module MOM_offline_aux

