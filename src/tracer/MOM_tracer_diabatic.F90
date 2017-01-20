!> This module contains routines that implement physical fluxes of tracers (e.g. due
!! to surface fluxes or mixing). These are intended to be called from call_tracer_column_fns
!! in the MOM_tracer_flow_control module.
module MOM_tracer_diabatic

! This file is part of MOM6. See LICENSE.md for the license.
use MOM_grid,             only : ocean_grid_type
use MOM_verticalGrid,     only : verticalGrid_type
use MOM_forcing_type,     only : forcing
use MOM_error_handler,    only : MOM_error, FATAL, WARNING

implicit none ; private

#include <MOM_memory.h>
public tracer_vertdiff
public applyTracerBoundaryFluxesInOut
!> This subroutine solves a tridiagonal equation for the final tracer
!! concentrations after the dual-entrainments, and possibly sinking or surface
!! and bottom sources, are applied.  The sinking is implemented with an
!! fully implicit upwind advection scheme.

contains

subroutine tracer_vertdiff(h_old, ea, eb, dt, tr, G, GV, &
                           sfc_flux, btm_flux, btm_reservoir, sink_rate, convert_flux_in)
  type(ocean_grid_type),                     intent(in)    :: G             !< ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV            !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_old         !< layer thickness before entrainment (m or kg m-2)
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: ea            !< amount of fluid entrained from the layer above (units of h_work)
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: eb            !< amount of fluid entrained from the layer below (units of h_work)
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: tr            !< tracer concentration (in concentration units CU)
  real,                                      intent(in)    :: dt            !< amount of time covered by this call (seconds)
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(in)    :: sfc_flux      !< surface flux of the tracer (in CU * kg m-2 s-1)
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(in)    :: btm_flux      !< The (negative upward) bottom flux of the tracer,
                                                                            !! in units of (CU * kg m-2 s-1)
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(inout) :: btm_reservoir !< amount of tracer in a bottom reservoir (units of CU kg m-2; formerly CU m)
  real,                             optional,intent(in)    :: sink_rate     !< rate at which the tracer sinks, in m s-1
  logical,                          optional,intent(in)  :: convert_flux_in    !< True if the specified sfc_flux needs to be integrated in time

  real :: sink_dist ! The distance the tracer sinks in a time step, in m or kg m-2.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    sfc_src, &      ! The time-integrated surface source of the tracer, in
                    ! units of m or kg m-2 times a concentration.
    btm_src         ! The time-integrated bottom source of the tracer, in
                    ! units of m or kg m-2  times a concentration.
  real, dimension(SZI_(G)) :: &
    b1, &           ! b1 is used by the tridiagonal solver, in m-1 or m2 kg-1.
    d1              ! d1=1-c1 is used by the tridiagonal solver, nondimensional.
  real :: c1(SZI_(G),SZK_(GV))    ! c1 is used by the tridiagonal solver, ND.
  real :: h_minus_dsink(SZI_(G),SZK_(GV)) ! The layer thickness minus the
                    ! difference in sinking rates across the layer, in m or kg m-2.
                    ! By construction, 0 <= h_minus_dsink < h_work.
  real :: sink(SZI_(G),SZK_(GV)+1) ! The tracer's sinking distances at the
                    ! interfaces, limited to prevent characteristics from
                    ! crossing within a single timestep, in m or kg m-2.
  real :: b_denom_1 ! The first term in the denominator of b1, in m or kg m-2.
  real :: h_tr      ! h_tr is h at tracer points with a h_neglect added to
                    ! ensure positive definiteness, in m or kg m-2.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in m.
  logical :: convert_flux = .true.


  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (present(convert_flux_in)) convert_flux = convert_flux_in
  h_neglect = GV%H_subroundoff
  sink_dist = 0.0
  if (present(sink_rate)) sink_dist = (dt*sink_rate) * GV%m_to_H
!$OMP parallel default(none) shared(is,ie,js,je,sfc_src,btm_src,sfc_flux,dt,G,GV,btm_flux, &
!$OMP                               sink_rate,btm_reservoir,nz,sink_dist,ea,      &
!$OMP                               h_old,convert_flux,h_neglect,eb,tr) &
!$OMP                       private(sink,h_minus_dsink,b_denom_1,b1,d1,h_tr,c1)
!$OMP do
  do j=js,je; do i=is,ie ; sfc_src(i,j) = 0.0 ; btm_src(i,j) = 0.0 ; enddo; enddo
  if (present(sfc_flux)) then
    if(convert_flux) then
!$OMP do
      do j = js, je; do i = is,ie
        sfc_src(i,j) = (sfc_flux(i,j)*dt) * GV%kg_m2_to_H
      enddo; enddo
    else
!$OMP do
      do j = js, je; do i = is,ie
        sfc_src(i,j) = sfc_flux(i,j)
      enddo; enddo
    endif
  endif
  if (present(btm_flux)) then
    if(convert_flux) then
!$OMP do
      do j = js, je; do i = is,ie
        btm_src(i,j) = (btm_flux(i,j)*dt) * GV%kg_m2_to_H
      enddo; enddo
    else
!$OMP do
      do j = js, je; do i = is,ie
        btm_src(i,j) = btm_flux(i,j)
      enddo; enddo
    endif
  endif

  if (present(sink_rate)) then
!$OMP do
    do j=js,je
      ! Find the sinking rates at all interfaces, limiting them if necesary
      ! so that the characteristics do not cross within a timestep.
      !   If a non-constant sinking rate were used, that would be incorprated
      ! here.
      if (present(btm_reservoir)) then
        do i=is,ie ; sink(i,nz+1) = sink_dist ; enddo
        do k=2,nz ; do i=is,ie
          sink(i,K) = sink_dist ; h_minus_dsink(i,k) = h_old(i,j,k)
        enddo ; enddo
      else
        do i=is,ie ; sink(i,nz+1) = 0.0 ; enddo
        ! Find the limited sinking distance at the interfaces.
        do k=nz,2,-1 ; do i=is,ie
          if (sink(i,K+1) >= sink_dist) then
            sink(i,K) = sink_dist
            h_minus_dsink(i,k) = h_old(i,j,k) + (sink(i,K+1) - sink(i,K))
          elseif (sink(i,K+1) + h_old(i,j,k) < sink_dist) then
            sink(i,K) = sink(i,K+1) + h_old(i,j,k)
            h_minus_dsink(i,k) = 0.0
          else
            sink(i,K) = sink_dist
            h_minus_dsink(i,k) = (h_old(i,j,k) + sink(i,K+1)) - sink(i,K)
          endif
        enddo ; enddo
      endif
      do i=is,ie
        sink(i,1) = 0.0 ; h_minus_dsink(i,1) = (h_old(i,j,1) + sink(i,2))
      enddo

      ! Now solve the tridiagonal equation for the tracer concentrations.
      do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
        b_denom_1 = h_minus_dsink(i,1) + ea(i,j,1) + h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,1))
        d1(i) = b_denom_1 * b1(i)
        h_tr = h_old(i,j,1) + h_neglect
        tr(i,j,1) = (b1(i)*h_tr)*tr(i,j,1) + sfc_src(i,j)
      endif ; enddo
      do k=2,nz-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
        c1(i,k) = eb(i,j,k-1) * b1(i)
        b_denom_1 = h_minus_dsink(i,k) + d1(i) * (ea(i,j,k) + sink(i,K)) + &
                    h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
        d1(i) = b_denom_1 * b1(i)
        h_tr = h_old(i,j,k) + h_neglect
        tr(i,j,k) = b1(i) * (h_tr * tr(i,j,k) + &
                             (ea(i,j,k) + sink(i,K)) * tr(i,j,k-1))
      endif ; enddo ; enddo
      do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
        c1(i,nz) = eb(i,j,nz-1) * b1(i)
        b_denom_1 = h_minus_dsink(i,nz) + d1(i) * (ea(i,j,nz) + sink(i,nz)) + &
                    h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,nz))
        h_tr = h_old(i,j,nz) + h_neglect
        tr(i,j,nz) = b1(i) * ((h_tr * tr(i,j,nz) + btm_src(i,j)) + &
                              (ea(i,j,nz) + sink(i,nz)) * tr(i,j,nz-1))
      endif ; enddo
      if (present(btm_reservoir)) then ; do i=is,ie ; if (G%mask2dT(i,j)>0.5) then
        btm_reservoir(i,j) = btm_reservoir(i,j) + &
                             (sink(i,nz+1)*tr(i,j,nz)) * GV%H_to_kg_m2
      endif ; enddo ; endif

      do k=nz-1,1,-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
        tr(i,j,k) = tr(i,j,k) + c1(i,k+1)*tr(i,j,k+1)
      endif ; enddo ; enddo
    enddo
  else
!$OMP do
    do j=js,je
      do i=is,ie ; if (G%mask2dT(i,j) > -0.5) then
        h_tr = h_old(i,j,1) + h_neglect
        b_denom_1 = h_tr + ea(i,j,1)
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,1))
        d1(i) = h_tr * b1(i)
        tr(i,j,1) = (b1(i)*h_tr)*tr(i,j,1) + sfc_src(i,j)
       endif
      enddo
      do k=2,nz-1 ; do i=is,ie ; if (G%mask2dT(i,j) > -0.5) then
        c1(i,k) = eb(i,j,k-1) * b1(i)
        h_tr = h_old(i,j,k) + h_neglect
        b_denom_1 = h_tr + d1(i) * ea(i,j,k)
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
        d1(i) = b_denom_1 * b1(i)
        tr(i,j,k) = b1(i) * (h_tr * tr(i,j,k) + ea(i,j,k) * tr(i,j,k-1))
      endif ; enddo ; enddo
      do i=is,ie ; if (G%mask2dT(i,j) > -0.5) then
        c1(i,nz) = eb(i,j,nz-1) * b1(i)
        h_tr = h_old(i,j,nz) + h_neglect
        b_denom_1 = h_tr + d1(i)*ea(i,j,nz)
        b1(i) = 1.0 / ( b_denom_1 + eb(i,j,nz))
        tr(i,j,nz) = b1(i) * (( h_tr * tr(i,j,nz) + btm_src(i,j)) + &
                              ea(i,j,nz) * tr(i,j,nz-1))
      endif ; enddo
      do k=nz-1,1,-1 ; do i=is,ie ; if (G%mask2dT(i,j) > -0.5) then
        tr(i,j,k) = tr(i,j,k) + c1(i,k+1)*tr(i,j,k+1)
      endif ; enddo ; enddo
    enddo
  endif

!$OMP end parallel

end subroutine tracer_vertdiff

subroutine applyTracerBoundaryFluxesInOut(G, GV, Tr, dt, fluxes, h, &
                                    evap_CFL_limit, minimum_forcing_depth, in_flux_optional, out_flux_optional)
! This routine is modeled after applyBoundaryFluxesInOut in MOM_diabatic_aux.F90
! NOTE: Please note that in this routine sfc_flux gets set to zero to ensure that the surface
!       flux of the tracer does not get applied again during a subsequent call to tracer_vertdif

  type(ocean_grid_type),                 intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),               intent(in)    :: GV        !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: Tr  !< Tracer concentration on T-cell
  real,                                  intent(in)    :: dt !< Time-step over which forcing is applied (s)
  type(forcing),                         intent(in) :: fluxes !< Surface fluxes container
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: h  !< Layer thickness in H units
  real,                                       intent(in)  :: evap_CFL_limit
  real,                                       intent(in)  :: minimum_forcing_depth
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: in_flux_optional ! The total time-integrated amount of tracer!
                                                                             ! that enters with freshwater
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: out_flux_optional ! The total time-integrated amount of tracer!
                                                                              ! that leaves with freshwater
  integer, parameter :: maxGroundings = 5
  integer :: numberOfGroundings, iGround(maxGroundings), jGround(maxGroundings)
  real :: H_limit_fluxes, IforcingDepthScale, Idt
  real :: dThickness, dTracer
  real :: fractionOfForcing, hOld, Ithickness
  real :: RivermixConst  ! A constant used in implementing river mixing, in Pa s.
  real, dimension(SZI_(G)) :: &
    netMassInOut, &  ! surface water fluxes (H units) over time step
    netMassIn,    &  ! mass entering ocean surface (H units) over a time step
    netMassOut       ! mass leaving ocean surface (H units) over a time step

  real, dimension(SZI_(G), SZK_(G))                     :: h2d, Tr2d
  real, dimension(SZI_(G),SZJ_(G))                      :: in_flux  ! The total time-integrated amount of tracer!
                                                                       ! that enters with freshwater
  real, dimension(SZI_(G),SZJ_(G))                      :: out_flux ! The total time-integrated amount of tracer!
                                                                        ! that leaves with freshwater
  real, dimension(SZI_(G))                              :: in_flux_1d, out_flux_1d
  real                                                  :: hGrounding(maxGroundings)
  real    :: Tr_in
  integer :: i, j, is, ie, js, je, k, nz, n, nsw
  character(len=45) :: mesg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

!  ! Only apply forcing if fluxes%sw is associated.
!  if (.not.ASSOCIATED(fluxes%sw)) return

  in_flux(:,:) = 0.0 ; out_flux(:,:) = 0.0
  if(present(in_flux_optional)) then
    do j=js,je ; do i=is,ie
      in_flux(i,j) = in_flux_optional(i,j)
    enddo; enddo
  endif
  if(present(out_flux_optional)) then
    do j=js,je ; do i=is,ie
      out_flux(i,j) = out_flux_optional(i,j)
    enddo ; enddo
  endif

  Idt = 1.0/dt
  numberOfGroundings = 0

!$OMP parallel do default(none) shared(is,ie,js,je,nz,h,Tr,G,GV,fluxes,dt,    &
!$OMP                                  IforcingDepthScale,minimum_forcing_depth, &
!$OMP                                  numberOfGroundings,iGround,jGround,      &
!$OMP                                  in_flux,out_flux,hGrounding,Idt,evap_CFL_limit) &
!$OMP                          private(h2d,Tr2d,netMassInOut,netMassOut,      &
!$OMP                                  in_flux_1d,out_flux_1d,fractionOfForcing,     &
!$OMP                                  dThickness,dTracer,hOld,Ithickness,           &
!$OMP                                  netMassIn, Tr_in)

  ! Work in vertical slices for efficiency
  do j=js,je

    ! Copy state into 2D-slice arrays
    do k=1,nz ; do i=is,ie
      h2d(i,k) = h(i,j,k)
      Tr2d(i,k) = Tr(i,j,k)
    enddo ; enddo

    do i = is,ie
      in_flux_1d(i) = in_flux(i,j)
      out_flux_1d(i) = out_flux(i,j)
    enddo
    ! The surface forcing is contained in the fluxes type.
    ! We aggregate the thermodynamic forcing for a time step into the following:
    ! These should have been set and stored during a call to applyBoundaryFluxesInOut
    ! netMassIn    = net mass entering at ocean surface over a timestep
    ! netMassOut   = net mass leaving ocean surface (H units) over a time step.
    !                netMassOut < 0 means mass leaves ocean.

    ! Note here that the aggregateFW flag has already been taken care of in the call to
    ! applyBoundaryFluxesInOut
    do i=is,ie
        netMassOut(i) = fluxes%netMassOut(i,j)
        netMassIn(i)  = fluxes%netMassIn(i,j)
    enddo

    ! Apply the surface boundary fluxes in three steps:
    ! A/ update concentration from mass entering the ocean
    ! B/ update concentration from mass leaving ocean.
    do i=is,ie
      if (G%mask2dT(i,j)>0.) then

        ! A/ Update tracer due to incoming mass flux.
        do k=1,1

          ! Change in state due to forcing
          dThickness = netMassIn(i) ! Since we are adding mass, we can use all of it
          dTracer = 0.

          ! Update the forcing by the part to be consumed within the present k-layer.
          ! If fractionOfForcing = 1, then updated netMassIn, netHeat, and netSalt vanish.
          netMassIn(i) = netMassIn(i) - dThickness
          dTracer = dTracer + in_flux_1d(i)
          in_flux_1d(i) = 0.0

          ! Update state
          hOld     = h2d(i,k)               ! Keep original thickness in hand
          h2d(i,k) = h2d(i,k) + dThickness  ! New thickness
          if (h2d(i,k) > 0.0) then
            Ithickness  = 1.0/h2d(i,k)      ! Inverse new thickness
            ! The "if"s below avoid changing T/S by roundoff unnecessarily
            if (dThickness /= 0. .or. dTracer /= 0.) tr2d(i,k) = (hOld*tr2d(i,k)+ dTracer)*Ithickness
          endif

        enddo ! k=1,1

        ! B/ Update tracer from mass leaving ocean
        do k=1,nz

          ! Place forcing into this layer if this layer has nontrivial thickness.
          ! For layers thin relative to 1/IforcingDepthScale, then distribute
          ! forcing into deeper layers.
          IforcingDepthScale = 1. / max(GV%H_subroundoff, minimum_forcing_depth*GV%m_to_H - netMassOut(i) )
          ! fractionOfForcing = 1.0, unless h2d is less than IforcingDepthScale.
          fractionOfForcing = min(1.0, h2d(i,k)*IforcingDepthScale)

          ! In the case with (-1)*netMassOut*fractionOfForcing greater than cfl*h, we
          ! limit the forcing applied to this cell, leaving the remaining forcing to
          ! be distributed downwards.
          if (-fractionOfForcing*netMassOut(i) > evap_CFL_limit*h2d(i,k)) then
            fractionOfForcing = -evap_CFL_limit*h2d(i,k)/netMassOut(i)
          endif

          ! Change in state due to forcing
          dThickness = max( fractionOfForcing*netMassOut(i), -h2d(i,k) )
          ! Note this is slightly different to how salt is currently treated
          dTracer =  fractionOfForcing*out_flux_1d(i)

          ! Update the forcing by the part to be consumed within the present k-layer.
          ! If fractionOfForcing = 1, then new netMassOut vanishes.
          netMassOut(i) = netMassOut(i) - dThickness
          out_flux_1d(i) = out_flux_1d(i) - dTracer

          ! Update state by the appropriate increment.
          hOld     = h2d(i,k)               ! Keep original thickness in hand
          h2d(i,k) = h2d(i,k) + dThickness  ! New thickness
          if (h2d(i,k) > 0.) then
            Ithickness  = 1.0/h2d(i,k) ! Inverse of new thickness
            Tr2d(i,k)    = (hOld*Tr2d(i,k) + dTracer)*Ithickness
          elseif (h2d(i,k) < 0.0) then ! h2d==0 is a special limit that needs no extra handling
            write(0,*) 'applyTracerBoundaryFluxesInOut(): lon,lat=',G%geoLonT(i,j),G%geoLatT(i,j)
            write(0,*) 'applyTracerBoundaryFluxesInOut(): netTr,netH=',in_flux_1d(i)-out_flux_1d(i),netMassInOut(i)
            write(0,*) 'applyTracerBoundaryFluxesInOut(): h(n),h(n+1),k=',hOld,h2d(i,k),k
            call MOM_error(FATAL, "MOM_tracer_vertical.F90, applyTracerBoundaryFluxesInOut(): "//&
                           "Complete mass loss in column!")
          endif

        enddo ! k

      ! Check if trying to apply fluxes over land points
      elseif((abs(in_flux_1d(i))+abs(out_flux_1d(i))+abs(netMassIn(i))+abs(netMassOut(i)))>0.) then
        write(0,*) 'applyTracerBoundaryFluxesInOut(): lon,lat=',G%geoLonT(i,j),G%geoLatT(i,j)
        write(0,*) 'applyTracerBoundaryFluxesInOut(): in_flux, out_flux, netMassIn,netMassOut=',&
                   in_flux_1d(i), out_flux_1d(i),netMassIn(i),netMassOut(i)
        call MOM_error(FATAL, "MOM_tracer_vertical.F90, applyTracerBoundaryFluxesInOut(): "//&
                              "Mass loss over land?")
      endif

      ! If anything remains after the k-loop, then we have grounded out, which is a problem.
      if (abs(in_flux_1d(i))+abs(out_flux_1d(i)) /= 0.0) then
!$OMP critical
        numberOfGroundings = numberOfGroundings +1
        if (numberOfGroundings<=maxGroundings) then
          iGround(numberOfGroundings) = i ! Record i,j location of event for
          jGround(numberOfGroundings) = j ! warning message
          hGrounding(numberOfGroundings) = abs(in_flux_1d(i))+abs(out_flux_1d(i))
        endif
!$OMP end critical
      endif

    enddo ! i

    ! Step C/ copy updated tracer concentration from the 2d slice now back into model state.
    do k=1,nz ; do i=is,ie
      Tr(i,j,k) = Tr2d(i,k)
      h(i,j,k) = h2d(i,k)
    enddo ; enddo

  enddo ! j-loop finish

  if (numberOfGroundings>0) then
    do i = 1, min(numberOfGroundings, maxGroundings)
      write(mesg(1:45),'(3es15.3)') G%geoLonT( iGround(i), jGround(i) ), &
                             G%geoLatT( iGround(i), jGround(i)) , hGrounding(i)
      call MOM_error(WARNING, "MOM_tracer_vertical.F90, applyTracerBoundaryFluxesInOut(): "//&
                              "Tracer created. x,y,dh= "//trim(mesg), all_print=.true.)
    enddo

    if (numberOfGroundings - maxGroundings > 0) then
      write(mesg, '(i4)') numberOfGroundings - maxGroundings
      call MOM_error(WARNING, "MOM_tracer_vertical.F90, applyTracerBoundaryFluxesInOut(): "//&
                              trim(mesg) // " groundings remaining")
    endif
  endif

end subroutine applyTracerBoundaryFluxesInOut
end module MOM_tracer_diabatic
