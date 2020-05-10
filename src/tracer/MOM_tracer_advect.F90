!>  This module contains the subroutines that advect tracers along coordinate surfaces.
module MOM_tracer_advect

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,       only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,       only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,   only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator,   only : register_diag_field, safe_alloc_ptr, time_type
use MOM_domains,         only : sum_across_PEs, max_across_PEs
use MOM_domains,         only : create_group_pass, do_group_pass, group_pass_type, pass_var
use MOM_error_handler,   only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,     only : get_param, log_version, param_file_type
use MOM_grid,            only : ocean_grid_type
use MOM_open_boundary,   only : ocean_OBC_type, OBC_NONE, OBC_DIRECTION_E
use MOM_open_boundary,   only : OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_open_boundary,   only : OBC_segment_type
use MOM_tracer_registry, only : tracer_registry_type, tracer_type
use MOM_unit_scaling,    only : unit_scale_type
use MOM_verticalGrid,    only : verticalGrid_type
implicit none ; private

#include <MOM_memory.h>

public advect_tracer
public tracer_advect_init
public tracer_advect_end

!> Control structure for this module
type, public :: tracer_advect_CS ; private
  real    :: dt                    !< The baroclinic dynamics time step [T ~> s].
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !< timing of diagnostic output.
  logical :: debug                 !< If true, write verbose checksums for debugging purposes.
  logical :: usePPM                !< If true, use PPM instead of PLM
  logical :: useHuynh              !< If true, use the Huynh scheme for PPM interface values
  type(group_pass_type) :: pass_uhr_vhr_t_hprev !< A structred used for group passes
end type tracer_advect_CS

!>@{ CPU time clocks
integer :: id_clock_advect
integer :: id_clock_pass
integer :: id_clock_sync
!>@}

contains

!> This routine time steps the tracer concentration using a
!! monotonic, conservative, weakly diffusive scheme.
subroutine advect_tracer(h_end, uhtr, vhtr, OBC, dt, G, GV, US, CS, Reg, &
      h_prev_opt, max_iter_in, x_first_in, uhr_out, vhr_out, h_out)
  type(ocean_grid_type),   intent(inout) :: G     !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV    !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: h_end !< layer thickness after advection [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: uhtr  !< accumulated volume/mass flux through zonal face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(in)    :: vhtr  !< accumulated volume/mass flux through merid face [H L2 ~> m3 or kg]
  type(ocean_OBC_type),    pointer       :: OBC   !< specifies whether, where, and what OBCs are used
  real,                    intent(in)    :: dt    !< time increment [T ~> s]
  type(unit_scale_type),   intent(in)    :: US    !< A dimensional unit scaling type
  type(tracer_advect_CS),  pointer       :: CS    !< control structure for module
  type(tracer_registry_type), pointer    :: Reg   !< pointer to tracer registry
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                 optional, intent(in)    :: h_prev_opt !< layer thickness before advection [H ~> m or kg m-2]
  integer,       optional, intent(in)    :: max_iter_in !< The maximum number of iterations
  logical,       optional, intent(in)    :: x_first_in !< If present, indicate whether to update
                                                  !! first in the x- or y-direction.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                 optional, intent(out)    :: uhr_out  !< accumulated volume/mass flux through zonal face
                                                  !! [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                 optional, intent(out)    :: vhr_out  !< accumulated volume/mass flux through merid face
                                                  !! [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                 optional, intent(out)    :: h_out !< layer thickness before advection [H ~> m or kg m-2]

  type(tracer_type) :: Tr(MAX_FIELDS_) ! The array of registered tracers
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    hprev           ! cell volume at the end of previous tracer change [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: &
    uhr             ! The remaining zonal thickness flux [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: &
    vhr             ! The remaining meridional thickness fluxes [H L2 ~> m3 or kg]
  real :: uh_neglect(SZIB_(G),SZJ_(G)) ! uh_neglect and vh_neglect are the
  real :: vh_neglect(SZI_(G),SZJB_(G)) ! magnitude of remaining transports that
                                       ! can be simply discarded [H L2 ~> m3 or kg].

  real :: landvolfill                   ! An arbitrary? nonzero cell volume [H L2 ~> m3 or kg].
  real :: Idt                           ! 1/dt [T-1 ~> s-1].
  logical :: domore_u(SZJ_(G),SZK_(G))  ! domore__ indicate whether there is more
  logical :: domore_v(SZJB_(G),SZK_(G)) ! advection to be done in the corresponding
                                        ! row or column.
  logical :: x_first            ! If true, advect in the x-direction first.
  integer :: max_iter           ! maximum number of iterations in each layer
  integer :: domore_k(SZK_(G))
  integer :: stencil            ! stencil of the advection scheme
  integer :: nsten_halo         ! number of stencils that fit in the halos
  integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz, itt, ntr, do_any
  integer :: isv, iev, jsv, jev ! The valid range of the indices.
  integer :: IsdB, IedB, JsdB, JedB

  domore_u(:,:) = .false.
  domore_v(:,:) = .false.
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  landvolfill = 1.0e-20         ! This is arbitrary, but must be positive.
  stencil = 2                   ! The scheme's stencil; 2 for PLM and PPM:H3

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_tracer_advect: "// &
       "tracer_advect_init must be called before advect_tracer.")
  if (.not. associated(Reg)) call MOM_error(FATAL, "MOM_tracer_advect: "// &
       "register_tracer must be called before advect_tracer.")
  if (Reg%ntr==0) return
  call cpu_clock_begin(id_clock_advect)
  x_first = (MOD(G%first_direction,2) == 0)

  ! increase stencil size for Colella & Woodward PPM
  if (CS%usePPM .and. .not. CS%useHuynh) stencil = 3

  ntr = Reg%ntr
  do m=1,ntr ; Tr(m) = Reg%Tr(m) ; enddo
  Idt = 1.0 / dt

  max_iter = 2*INT(CEILING(dt/CS%dt)) + 1

  if (present(max_iter_in)) max_iter = max_iter_in
  if (present(x_first_in))  x_first = x_first_in
  call cpu_clock_begin(id_clock_pass)
  call create_group_pass(CS%pass_uhr_vhr_t_hprev, uhr, vhr, G%Domain)
  call create_group_pass(CS%pass_uhr_vhr_t_hprev, hprev, G%Domain)
  do m=1,ntr
    call create_group_pass(CS%pass_uhr_vhr_t_hprev, Tr(m)%t, G%Domain)
  enddo
  call cpu_clock_end(id_clock_pass)

!$OMP parallel default(none) shared(nz,jsd,jed,IsdB,IedB,uhr,jsdB,jedB,Isd,Ied,vhr, &
!$OMP                               hprev,domore_k,js,je,is,ie,uhtr,vhtr,G,GV,h_end,&
!$OMP                               uh_neglect,vh_neglect,ntr,Tr,h_prev_opt)

! This initializes the halos of uhr and vhr because pass_vector might do
! calculations on them, even though they are never used.
!$OMP do

  do k=1,nz
    do j=jsd,jed ; do I=IsdB,IedB ; uhr(I,j,k) = 0.0 ; enddo ; enddo
    do J=jsdB,jedB ; do i=Isd,Ied ; vhr(i,J,k) = 0.0 ; enddo ; enddo
    do j=jsd,jed ; do i=Isd,Ied ; hprev(i,j,k) = 0.0 ; enddo ; enddo
    domore_k(k)=1
!  Put the remaining (total) thickness fluxes into uhr and vhr.
    do j=js,je ; do I=is-1,ie ; uhr(I,j,k) = uhtr(I,j,k) ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; vhr(i,J,k) = vhtr(i,J,k) ; enddo ; enddo
    if (.not. present(h_prev_opt)) then
    !   This loop reconstructs the thickness field the last time that the
    ! tracers were updated, probably just after the diabatic forcing.  A useful
    ! diagnostic could be to compare this reconstruction with that older value.
      do i=is,ie ; do j=js,je
        hprev(i,j,k) = max(0.0, G%areaT(i,j)*h_end(i,j,k) + &
             ((uhr(I,j,k) - uhr(I-1,j,k)) + (vhr(i,J,k) - vhr(i,J-1,k))))
    ! In the case that the layer is now dramatically thinner than it was previously,
    ! add a bit of mass to avoid truncation errors.  This will lead to
    ! non-conservation of tracers
        hprev(i,j,k) = hprev(i,j,k) + &
                       max(0.0, 1.0e-13*hprev(i,j,k) - G%areaT(i,j)*h_end(i,j,k))
      enddo ; enddo
    else
      do i=is,ie ; do j=js,je
        hprev(i,j,k) = h_prev_opt(i,j,k)
      enddo ; enddo
    endif
  enddo


!$OMP do
  do j=jsd,jed ; do I=isd,ied-1
    uh_neglect(I,j) = GV%H_subroundoff*MIN(G%areaT(i,j),G%areaT(i+1,j))
  enddo ; enddo
!$OMP do
  do J=jsd,jed-1 ; do i=isd,ied
    vh_neglect(i,J) = GV%H_subroundoff*MIN(G%areaT(i,j),G%areaT(i,j+1))
  enddo ; enddo

!$OMP do
  ! initialize diagnostic fluxes and tendencies
  do m=1,ntr
    if (associated(Tr(m)%ad_x)) then
      do k=1,nz ; do j=jsd,jed ; do i=isd,ied
        Tr(m)%ad_x(I,j,k) = 0.0
      enddo ; enddo ; enddo
    endif
    if (associated(Tr(m)%ad_y)) then
      do k=1,nz ; do J=jsd,jed ; do i=isd,ied
        Tr(m)%ad_y(i,J,k) = 0.0
      enddo ; enddo ; enddo
    endif
    if (associated(Tr(m)%advection_xy)) then
      do k=1,nz ; do j=jsd,jed ; do i=isd,ied
        Tr(m)%advection_xy(i,j,k) = 0.0
      enddo ; enddo ; enddo
    endif
    if (associated(Tr(m)%ad2d_x)) then
      do j=jsd,jed ; do i=isd,ied ; Tr(m)%ad2d_x(I,j) = 0.0 ; enddo ; enddo
    endif
    if (associated(Tr(m)%ad2d_y)) then
      do J=jsd,jed ; do i=isd,ied ; Tr(m)%ad2d_y(i,J) = 0.0 ; enddo ; enddo
    endif
  enddo
!$OMP end parallel

  isv = is ; iev = ie ; jsv = js ; jev = je

  do itt=1,max_iter

    if (isv > is-stencil) then
      call do_group_pass(CS%pass_uhr_vhr_t_hprev, G%Domain, clock=id_clock_pass)

      nsten_halo = min(is-isd,ied-ie,js-jsd,jed-je)/stencil
      isv = is-nsten_halo*stencil ; jsv = js-nsten_halo*stencil
      iev = ie+nsten_halo*stencil ; jev = je+nsten_halo*stencil
      ! Reevaluate domore_u & domore_v unless the valid range is the same size as
      ! before.  Also, do this if there is Strang splitting.
      if ((nsten_halo > 1) .or. (itt==1)) then
!$OMP parallel do default(none) shared(nz,domore_k,jsv,jev,domore_u,isv,iev,stencil, &
!$OMP                                  uhr,domore_v,vhr)
        do k=1,nz ; if (domore_k(k) > 0) then
          do j=jsv,jev ; if (.not.domore_u(j,k)) then
            do i=isv+stencil-1,iev-stencil; if (uhr(I,j,k) /= 0.0) then
              domore_u(j,k) = .true. ; exit
            endif ; enddo ! i-loop
          endif ; enddo
          do J=jsv+stencil-1,jev-stencil ; if (.not.domore_v(J,k)) then
            do i=isv+stencil,iev-stencil; if (vhr(i,J,k) /= 0.0) then
              domore_v(J,k) = .true. ; exit
            endif ; enddo ! i-loop
          endif ; enddo

          !   At this point, domore_k is global.  Change it so that it indicates
          ! whether any work is needed on a layer on this processor.
          domore_k(k) = 0
          do j=jsv,jev ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
          do J=jsv+stencil-1,jev-stencil ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo

        endif ; enddo ! k-loop
      endif
    endif

    ! Set the range of valid points after this iteration.
    isv = isv + stencil ; iev = iev - stencil
    jsv = jsv + stencil ; jev = jev - stencil

    !  To ensure positive definiteness of the thickness at each iteration, the
    !  mass fluxes out of each layer are checked each step, and limited to keep
    !  the thicknesses positive.  This means that several iterations may be required
    !  for all the transport to happen.  The sum over domore_k keeps the processors
    !  synchronized.  This may not be very efficient, but it should be reliable.

!$OMP parallel default(private) shared(nz,domore_k,x_first,Tr,hprev,uhr,uh_neglect,  &
!$OMP                                  OBC,domore_u,ntr,Idt,isv,iev,jsv,jev,stencil, &
!$OMP                                  G,GV,CS,vhr,vh_neglect,domore_v,US)

    if (x_first) then

      !$OMP do ordered
      do k=1,nz ; if (domore_k(k) > 0) then
        ! First, advect zonally.
        call advect_x(Tr, hprev, uhr, uh_neglect, OBC, domore_u, ntr, Idt, &
                      isv, iev, jsv-stencil, jev+stencil, k, G, GV, US, CS%usePPM, CS%useHuynh)
      endif ; enddo

      !$OMP do ordered
      do k=1,nz ; if (domore_k(k) > 0) then
        !  Next, advect meridionally.
        call advect_y(Tr, hprev, vhr, vh_neglect, OBC, domore_v, ntr, Idt, &
                      isv, iev, jsv, jev, k, G, GV, US, CS%usePPM, CS%useHuynh)

        ! Update domore_k(k) for the next iteration
        domore_k(k) = 0
        do j=jsv-stencil,jev+stencil ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
        do J=jsv-1,jev ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo

      endif ; enddo

    else

      !$OMP do ordered
      do k=1,nz ; if (domore_k(k) > 0) then
        ! First, advect meridionally.
        call advect_y(Tr, hprev, vhr, vh_neglect, OBC, domore_v, ntr, Idt, &
                      isv-stencil, iev+stencil, jsv, jev, k, G, GV, US, CS%usePPM, CS%useHuynh)
      endif ; enddo

      !$OMP do ordered
      do k=1,nz ; if (domore_k(k) > 0) then
        ! Next, advect zonally.
        call advect_x(Tr, hprev, uhr, uh_neglect, OBC, domore_u, ntr, Idt, &
                      isv, iev, jsv, jev, k, G, GV, US, CS%usePPM, CS%useHuynh)

        ! Update domore_k(k) for the next iteration
        domore_k(k) = 0
        do j=jsv,jev ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
        do J=jsv-1,jev ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo
      endif ; enddo

    endif ! x_first

!$OMP end parallel

    ! If the advection just isn't finishing after max_iter, move on.
    if (itt >= max_iter) then
      exit
    endif

    ! Exit if there are no layers that need more iterations.
    if (isv > is-stencil) then
      do_any = 0
      call cpu_clock_begin(id_clock_sync)
      call sum_across_PEs(domore_k(:), nz)
      call cpu_clock_end(id_clock_sync)
      do k=1,nz ; do_any = do_any + domore_k(k) ; enddo
      if (do_any == 0) then
        exit
      endif

    endif

  enddo ! Iterations loop

  if (present(uhr_out)) uhr_out(:,:,:) = uhr(:,:,:)
  if (present(vhr_out)) vhr_out(:,:,:) = vhr(:,:,:)
  if (present(h_out)) h_out(:,:,:) = hprev(:,:,:)

  call cpu_clock_end(id_clock_advect)

end subroutine advect_tracer


!> This subroutine does 1-d flux-form advection in the zonal direction using
!! a monotonic piecewise linear scheme.
subroutine advect_x(Tr, hprev, uhr, uh_neglect, OBC, domore_u, ntr, Idt, &
                    is, ie, js, je, k, G, GV, US, usePPM, useHuynh)
  type(ocean_grid_type),                     intent(inout) :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)    :: GV   !< The ocean's vertical grid structure
  type(tracer_type), dimension(ntr),         intent(inout) :: Tr   !< The array of registered tracers to work on
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: hprev !< cell volume at the end of previous
                                                                  !! tracer change [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uhr !< accumulated volume/mass flux through
                                                                  !! the zonal face [H L2 ~> m3 or kg]
   real, dimension(SZIB_(G),SZJ_(G)),        intent(in)    :: uh_neglect !< A tiny zonal mass flux that can
                                                                  !! be neglected [H L2 ~> m3 or kg]
  type(ocean_OBC_type),                      pointer       :: OBC !< specifies whether, where, and what OBCs are used
  logical, dimension(SZJ_(G),SZK_(G)),       intent(inout) :: domore_u !< If true, there is more advection to be
                                                                  !! done in this u-row
  real,                                      intent(in)    :: Idt !< The inverse of dt [T-1 ~> s-1]
  integer,                                   intent(in)    :: ntr !< The number of tracers
  integer,                                   intent(in)    :: is  !< The starting tracer i-index to work on
  integer,                                   intent(in)    :: ie  !< The ending tracer i-index to work on
  integer,                                   intent(in)    :: js  !< The starting tracer j-index to work on
  integer,                                   intent(in)    :: je  !< The ending tracer j-index to work on
  integer,                                   intent(in)    :: k   !< The k-level to work on
  type(unit_scale_type),                     intent(in)    :: US  !< A dimensional unit scaling type
  logical,                                   intent(in)    :: usePPM !< If true, use PPM instead of PLM
  logical,                                   intent(in)    :: useHuynh !< If true, use the Huynh scheme
                                                                     !! for PPM interface values

  real, dimension(SZI_(G),ntr) :: &
    slope_x             ! The concentration slope per grid point [conc].
  real, dimension(SZIB_(G),SZJ_(G),ntr) :: &
    flux_x              ! The tracer flux across a boundary [H L2 conc ~> m3 conc or kg conc].
  real, dimension(SZI_(G),ntr) :: &
    T_tmp               ! The copy of the tracer concentration at constant i,k [H m2 conc ~> m3 conc or kg conc].

  real :: maxslope      ! The maximum concentration slope per grid point
                        ! consistent with monotonicity [conc].
  real :: hup, hlos     ! hup is the upwind volume, hlos is the
                        ! part of that volume that might be lost
                        ! due to advection out the other side of
                        ! the grid box, both in [H L2 ~> m3 or kg].
  real :: uhh(SZIB_(G)) ! The zonal flux that occurs during the
                        ! current iteration [H L2 ~> m3 or kg].
  real, dimension(SZIB_(G)) :: &
    hlst, &             ! Work variable [H L2 ~> m3 or kg].
    Ihnew, &            ! Work variable [H-1 L-2 ~> m-3 or kg-1].
    CFL                 ! A nondimensional work variable [nondim].
  real :: min_h         ! The minimum thickness that can be realized during
                        ! any of the passes [H ~> m or kg m-2].
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
  logical :: do_i(SZIB_(G),SZJ_(G))     ! If true, work on given points.
  logical :: do_any_i
  integer :: i, j, m, n, i_up, stencil
  real :: aR, aL, dMx, dMn, Tp, Tc, Tm, dA, mA, a6
  real :: fac1,u_L_in,u_L_out  ! terms used for time-stepping OBC reservoirs
  type(OBC_segment_type), pointer :: segment=>NULL()
  logical :: usePLMslope
  logical, dimension(SZJ_(G),SZK_(G)) :: domore_u_initial

  ! keep a local copy of the initial values of domore_u, which is to be used when computing ad2d_x
  ! diagnostic at the end of this subroutine.
  domore_u_initial = domore_u

  usePLMslope = .not. (usePPM .and. useHuynh)
  ! stencil for calculating slope values
  stencil = 1
  if (usePPM .and. .not. useHuynh) stencil = 2

  min_h = 0.1*GV%Angstrom_H
  h_neglect = GV%H_subroundoff

! do I=is-1,ie ; ts2(I) = 0.0 ; enddo
  do I=is-1,ie ; CFL(I) = 0.0 ; enddo

  do j=js,je ; if (domore_u(j,k)) then
    domore_u(j,k) = .false.

    ! Calculate the i-direction profiles (slopes) of each tracer that
    ! is being advected.
    if (usePLMslope) then
      do m=1,ntr ; do i=is-stencil,ie+stencil
       !if (ABS(Tr(m)%t(i+1,j,k)-Tr(m)%t(i,j,k)) < &
       !    ABS(Tr(m)%t(i,j,k)-Tr(m)%t(i-1,j,k))) then
       !  maxslope = 4.0*(Tr(m)%t(i+1,j,k)-Tr(m)%t(i,j,k))
       !else
       !  maxslope = 4.0*(Tr(m)%t(i,j,k)-Tr(m)%t(i-1,j,k))
       !endif
       !if ((Tr(m)%t(i+1,j,k)-Tr(m)%t(i,j,k)) * (Tr(m)%t(i,j,k)-Tr(m)%t(i-1,j,k)) < 0.0) then
       !  slope_x(i,m) = 0.0
       !elseif (ABS(Tr(m)%t(i+1,j,k)-Tr(m)%t(i-1,j,k))<ABS(maxslope)) then
       !  slope_x(i,m) = G%mask2dCu(I,j)*G%mask2dCu(I-1,j) * &
       !                 0.5*(Tr(m)%t(i+1,j,k)-Tr(m)%t(i-1,j,k))
       !else
       !  slope_x(i,m) = G%mask2dCu(I,j)*G%mask2dCu(I-1,j) * 0.5*maxslope
       !endif
        Tp = Tr(m)%t(i+1,j,k) ; Tc = Tr(m)%t(i,j,k) ; Tm = Tr(m)%t(i-1,j,k)
        dMx = max( Tp, Tc, Tm ) - Tc
        dMn= Tc - min( Tp, Tc, Tm )
        slope_x(i,m) = G%mask2dCu(I,j)*G%mask2dCu(I-1,j) * &
            sign( min(0.5*abs(Tp-Tm), 2.0*dMx, 2.0*dMn), Tp-Tm )
      enddo ; enddo
    endif ! usePLMslope

    ! make a copy of the tracers in case values need to be overridden for OBCs
    do m = 1,ntr
      do i=G%isd,G%ied
        T_tmp(i,m) = Tr(m)%t(i,j,k)
      enddo
    enddo
    ! loop through open boundaries and recalculate flux terms
    if (associated(OBC)) then ; if (OBC%OBC_pe) then
       do n=1,OBC%number_of_segments
         segment=>OBC%segment(n)
         if (.not. associated(segment%tr_Reg)) cycle
         if (segment%is_E_or_W) then
           if (j>=segment%HI%jsd .and. j<=segment%HI%jed) then
              I = segment%HI%IsdB
              do m = 1,ntr ! replace tracers with OBC values
                if (associated(segment%tr_Reg%Tr(m)%tres)) then
                   if (segment%direction == OBC_DIRECTION_W) then
                      T_tmp(i,m) = segment%tr_Reg%Tr(m)%tres(i,j,k)
                   else
                      T_tmp(I+1,m) = segment%tr_Reg%Tr(m)%tres(i,j,k)
                   endif
                else
                   if (segment%direction == OBC_DIRECTION_W) then
                      T_tmp(i,m) = segment%tr_Reg%Tr(m)%OBC_inflow_conc
                   else
                      T_tmp(I+1,m) = segment%tr_Reg%Tr(m)%OBC_inflow_conc
                   endif
                endif
              enddo
              do m = 1,ntr ! Apply update tracer values for slope calculation
                do i=segment%HI%IsdB-1,segment%HI%IsdB+1
                  Tp = T_tmp(i+1,m) ; Tc = T_tmp(i,m) ; Tm = T_tmp(i-1,m)
                  dMx = max( Tp, Tc, Tm ) - Tc
                  dMn= Tc - min( Tp, Tc, Tm )
                  slope_x(i,m) = G%mask2dCu(I,j)*G%mask2dCu(I-1,j) * &
                       sign( min(0.5*abs(Tp-Tm), 2.0*dMx, 2.0*dMn), Tp-Tm )
                enddo
              enddo

           endif
         endif
       enddo
    endif; endif


    ! Calculate the i-direction fluxes of each tracer, using as much
    ! the minimum of the remaining mass flux (uhr) and the half the mass
    ! in the cell plus whatever part of its half of the mass flux that
    ! the flux through the other side does not require.
    do I=is-1,ie
      if (uhr(I,j,k) == 0.0) then
        uhh(I) = 0.0
        CFL(I) = 0.0
      elseif (uhr(I,j,k) < 0.0) then
        hup = hprev(i+1,j,k) - G%areaT(i+1,j)*min_h
        hlos = MAX(0.0,uhr(I+1,j,k))
        if ((((hup - hlos) + uhr(I,j,k)) < 0.0) .and. &
            ((0.5*hup + uhr(I,j,k)) < 0.0)) then
          uhh(I) = MIN(-0.5*hup,-hup+hlos,0.0)
          domore_u(j,k) = .true.
        else
          uhh(I) = uhr(I,j,k)
        endif
       !ts2(I) = 0.5*(1.0 + uhh(I) / (hprev(i+1,j,k) + h_neglect*G%areaT(i+1,j)))
        CFL(I) = - uhh(I) / (hprev(i+1,j,k) + h_neglect*G%areaT(i+1,j)) ! CFL is positive
      else
        hup = hprev(i,j,k) - G%areaT(i,j)*min_h
        hlos = MAX(0.0,-uhr(I-1,j,k))
        if ((((hup - hlos) - uhr(I,j,k)) < 0.0) .and. &
            ((0.5*hup - uhr(I,j,k)) < 0.0)) then
          uhh(I) = MAX(0.5*hup,hup-hlos,0.0)
          domore_u(j,k) = .true.
        else
          uhh(I) = uhr(I,j,k)
        endif
       !ts2(I) = 0.5*(1.0 - uhh(I) / (hprev(i,j,k) + h_neglect*G%areaT(i,j)))
        CFL(I) = uhh(I) / (hprev(i,j,k) + h_neglect*G%areaT(i,j)) ! CFL is positive
      endif
    enddo


    if (usePPM) then
      do m=1,ntr ; do I=is-1,ie
        ! centre cell depending on upstream direction
        if (uhh(I) >= 0.0) then
          i_up = i
        else
          i_up = i+1
        endif

        ! Implementation of PPM-H3
        Tp = T_tmp(i_up+1,m) ; Tc = T_tmp(i_up,m) ; Tm = T_tmp(i_up-1,m)

        if (useHuynh) then
          aL = ( 5.*Tc + ( 2.*Tm - Tp ) )/6. ! H3 estimate
          aL = max( min(Tc,Tm), aL) ; aL = min( max(Tc,Tm), aL) ! Bound
          aR = ( 5.*Tc + ( 2.*Tp - Tm ) )/6. ! H3 estimate
          aR = max( min(Tc,Tp), aR) ; aR = min( max(Tc,Tp), aR) ! Bound
        else
          aL = 0.5 * ((Tm + Tc) + (slope_x(i_up-1,m) - slope_x(i_up,m)) / 3.)
          aR = 0.5 * ((Tc + Tp) + (slope_x(i_up,m) - slope_x(i_up+1,m)) / 3.)
        endif

        dA = aR - aL ; mA = 0.5*( aR + aL )
        if (G%mask2dCu(I_up,j)*G%mask2dCu(I_up-1,j)*(Tp-Tc)*(Tc-Tm) <= 0.) then
           aL = Tc ; aR = Tc ! PCM for local extremum and bounadry cells
        elseif ( dA*(Tc-mA) > (dA*dA)/6. ) then
           aL = 3.*Tc - 2.*aR
        elseif ( dA*(Tc-mA) < - (dA*dA)/6. ) then
           aR = 3.*Tc - 2.*aL
        endif

        a6 = 6.*Tc - 3. * (aR + aL) ! Curvature

        if (uhh(I) >= 0.0) then
          flux_x(I,j,m) = uhh(I)*( aR - 0.5 * CFL(I) * ( &
               ( aR - aL ) - a6 * ( 1. - 2./3. * CFL(I) ) ) )
        else
          flux_x(I,j,m) = uhh(I)*( aL + 0.5 * CFL(I) * ( &
               ( aR - aL ) + a6 * ( 1. - 2./3. * CFL(I) ) ) )
        endif
      enddo ; enddo
    else ! PLM
      do m=1,ntr ; do I=is-1,ie
        if (uhh(I) >= 0.0) then
          ! Indirect implementation of PLM
         !aL = Tr(m)%t(i,j,k) - 0.5 * slope_x(i,m)
         !aR = Tr(m)%t(i,j,k) + 0.5 * slope_x(i,m)
         !flux_x(I,j,m) = uhh(I)*( aR - 0.5 * (aR-aL) * CFL(I) )
          ! Alternative implementation of PLM
         !aR = Tr(m)%t(i,j,k) + 0.5 * slope_x(i,m)
         !flux_x(I,j,m) = uhh(I)*( aR - 0.5 * slope_x(i,m) * CFL(I) )
          ! Alternative implementation of PLM
          Tc = T_tmp(i,m)
          flux_x(I,j,m) = uhh(I)*( Tc + 0.5 * slope_x(i,m) * ( 1. - CFL(I) ) )
          ! Original implementation of PLM
         !flux_x(I,j,m) = uhh(I)*(Tr(m)%t(i,j,k) + slope_x(i,m)*ts2(I))
        else
          ! Indirect implementation of PLM
         !aL = Tr(m)%t(i+1,j,k) - 0.5 * slope_x(i+1,m)
         !aR = Tr(m)%t(i+1,j,k) + 0.5 * slope_x(i+1,m)
         !flux_x(I,j,m) = uhh(I)*( aL + 0.5 * (aR-aL) * CFL(I) )
          ! Alternative implementation of PLM
         !aL = Tr(m)%t(i+1,j,k) - 0.5 * slope_x(i+1,m)
         !flux_x(I,j,m) = uhh(I)*( aL + 0.5 * slope_x(i+1,m) * CFL(I) )
          ! Alternative implementation of PLM
          Tc = T_tmp(i+1,m)
          flux_x(I,j,m) = uhh(I)*( Tc - 0.5 * slope_x(i+1,m) * ( 1. - CFL(I) ) )
          ! Original implementation of PLM
         !flux_x(I,j,m) = uhh(I)*(Tr(m)%t(i+1,j,k) - slope_x(i+1,m)*ts2(I))
        endif
       !ts2(I) = 0.5*(1.0 - uhh(I)/(hprev(i,j,k)+h_neglect*G%areaT(i,j)))
      enddo ; enddo
    endif ! usePPM

    if (associated(OBC)) then ; if (OBC%OBC_pe) then
      if (OBC%specified_u_BCs_exist_globally .or. OBC%open_u_BCs_exist_globally) then
        do n=1,OBC%number_of_segments
          segment=>OBC%segment(n)
          if (.not. associated(segment%tr_Reg)) cycle
          if (segment%is_E_or_W) then
            if (j>=segment%HI%jsd .and. j<=segment%HI%jed) then
              I = segment%HI%IsdB
              ! Tracer fluxes are set to prescribed values only for inflows from masked areas.
              ! Now changing to simply fixed inflows.
              if ((uhr(I,j,k) > 0.0) .and. (segment%direction == OBC_DIRECTION_W) .or. &
                 (uhr(I,j,k) < 0.0) .and. (segment%direction == OBC_DIRECTION_E)) then
                uhh(I) = uhr(I,j,k)
              ! should the reservoir evolve for this case Kate ?? - Nope
                do m=1,ntr
                  if (associated(segment%tr_Reg%Tr(m)%tres)) then
                    flux_x(I,j,m) = uhh(I)*segment%tr_Reg%Tr(m)%tres(I,j,k)
                  else ; flux_x(I,j,m) = uhh(I)*segment%tr_Reg%Tr(m)%OBC_inflow_conc ; endif
                enddo
              endif
            endif
          endif
        enddo
      endif

      if (OBC%open_u_BCs_exist_globally) then
        do n=1,OBC%number_of_segments
          segment=>OBC%segment(n)
          I = segment%HI%IsdB
          if (segment%is_E_or_W .and. (j >= segment%HI%jsd .and. j<= segment%HI%jed)) then
            if (segment%specified) cycle
            if (.not. associated(segment%tr_Reg)) cycle

            ! Tracer fluxes are set to prescribed values only for inflows from masked areas.
            if ((uhr(I,j,k) > 0.0) .and. (G%mask2dT(i,j) < 0.5) .or. &
               (uhr(I,j,k) < 0.0) .and. (G%mask2dT(i+1,j) < 0.5)) then
              uhh(I) = uhr(I,j,k)
              do m=1,ntr
                if (associated(segment%tr_Reg%Tr(m)%tres)) then
                  flux_x(I,j,m) = uhh(I)*segment%tr_Reg%Tr(m)%tres(I,j,k)
                else; flux_x(I,j,m) = uhh(I)*segment%tr_Reg%Tr(m)%OBC_inflow_conc; endif
              enddo
            endif
          endif
        enddo
      endif
    endif ; endif

    ! Calculate new tracer concentration in each cell after accounting
    ! for the i-direction fluxes.
    do I=is-1,ie
      uhr(I,j,k) = uhr(I,j,k) - uhh(I)
      if (abs(uhr(I,j,k)) < uh_neglect(I,j)) uhr(I,j,k) = 0.0
    enddo
    do i=is,ie
      if ((uhh(I) /= 0.0) .or. (uhh(I-1) /= 0.0)) then
        do_i(i,j) = .true.
        hlst(i) = hprev(i,j,k)
        hprev(i,j,k) = hprev(i,j,k) - (uhh(I) - uhh(I-1))
        if (hprev(i,j,k) <= 0.0) then ; do_i(i,j) = .false.
        elseif (hprev(i,j,k) < h_neglect*G%areaT(i,j)) then
          hlst(i) = hlst(i) + (h_neglect*G%areaT(i,j) - hprev(i,j,k))
          Ihnew(i) = 1.0 / (h_neglect*G%areaT(i,j))
        else ;  Ihnew(i) = 1.0 / hprev(i,j,k) ; endif
      else
        do_i(i,j) = .false.
      endif
    enddo

    ! update tracer concentration from i-flux and save some diagnostics
    do m=1,ntr

      ! update tracer
      do i=is,ie
        if (do_i(i,j)) then
          if (Ihnew(i) > 0.0) then
            Tr(m)%t(i,j,k) = (Tr(m)%t(i,j,k) * hlst(i) - &
                              (flux_x(I,j,m) - flux_x(I-1,j,m))) * Ihnew(i)
          endif
        endif
      enddo

      ! diagnostics
      if (associated(Tr(m)%ad_x)) then ; do i=is,ie ; if (do_i(i,j)) then
        Tr(m)%ad_x(I,j,k) = Tr(m)%ad_x(I,j,k) + flux_x(I,j,m)*Idt
      endif ; enddo ; endif

      ! diagnose convergence of flux_x (do not use the Ihnew(i) part of the logic).
      ! division by areaT to get into W/m2 for heat and kg/(s*m2) for salt.
      if (associated(Tr(m)%advection_xy)) then
        do i=is,ie ; if (do_i(i,j)) then
          Tr(m)%advection_xy(i,j,k) = Tr(m)%advection_xy(i,j,k) - (flux_x(I,j,m) - flux_x(I-1,j,m)) * &
                                          Idt * G%IareaT(i,j)
        endif ; enddo
      endif

    enddo

  endif


  enddo ! End of j-loop.

  ! compute ad2d_x diagnostic outside above j-loop so as to make the summation ordered when OMP is active.

  !$OMP ordered
  do j=js,je ; if (domore_u_initial(j,k)) then
    do m=1,ntr
      if (associated(Tr(m)%ad2d_x)) then ; do i=is,ie ; if (do_i(i,j)) then
        Tr(m)%ad2d_x(I,j) = Tr(m)%ad2d_x(I,j) + flux_x(I,j,m)*Idt
      endif ; enddo ; endif
    enddo
  endif ; enddo ! End of j-loop.
  !$OMP end ordered

end subroutine advect_x

!> This subroutine does 1-d flux-form advection using a monotonic piecewise
!! linear scheme.
subroutine advect_y(Tr, hprev, vhr, vh_neglect, OBC, domore_v, ntr, Idt, &
                    is, ie, js, je, k, G, GV, US, usePPM, useHuynh)
  type(ocean_grid_type),                     intent(inout) :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)    :: GV   !< The ocean's vertical grid structure
  type(tracer_type), dimension(ntr),         intent(inout) :: Tr   !< The array of registered tracers to work on
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: hprev !< cell volume at the end of previous
                                                                  !! tracer change [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vhr !< accumulated volume/mass flux through
                                                                  !! the meridional face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G)),         intent(inout) :: vh_neglect !< A tiny meridional mass flux that can
                                                                  !! be neglected [H L2 ~> m3 or kg]
  type(ocean_OBC_type),                      pointer       :: OBC !< specifies whether, where, and what OBCs are used
  logical, dimension(SZJB_(G),SZK_(G)),      intent(inout) :: domore_v !< If true, there is more advection to be
                                                                  !! done in this v-row
  real,                                      intent(in)    :: Idt !< The inverse of dt [T-1 ~> s-1]
  integer,                                   intent(in)    :: ntr !< The number of tracers
  integer,                                   intent(in)    :: is  !< The starting tracer i-index to work on
  integer,                                   intent(in)    :: ie  !< The ending tracer i-index to work on
  integer,                                   intent(in)    :: js  !< The starting tracer j-index to work on
  integer,                                   intent(in)    :: je  !< The ending tracer j-index to work on
  integer,                                   intent(in)    :: k   !< The k-level to work on
  type(unit_scale_type),                     intent(in)    :: US  !< A dimensional unit scaling type
  logical,                                   intent(in)    :: usePPM !< If true, use PPM instead of PLM
  logical,                                   intent(in)    :: useHuynh !< If true, use the Huynh scheme
                                                                     !! for PPM interface values

  real, dimension(SZI_(G),ntr,SZJ_(G)) :: &
    slope_y                     ! The concentration slope per grid point [conc].
  real, dimension(SZI_(G),ntr,SZJB_(G)) :: &
       flux_y                      ! The tracer flux across a boundary [H m2 conc ~> m3 conc or kg conc].
  real, dimension(SZI_(G),ntr,SZJB_(G)) :: &
    T_tmp               ! The copy of the tracer concentration at constant i,k [H m2 conc ~> m3 conc or kg conc].
  real :: maxslope              ! The maximum concentration slope per grid point
                                ! consistent with monotonicity [conc].
  real :: vhh(SZI_(G),SZJB_(G)) ! The meridional flux that occurs during the
                                ! current iteration [H L2 ~> m3 or kg].
  real :: hup, hlos             ! hup is the upwind volume, hlos is the
                                ! part of that volume that might be lost
                                ! due to advection out the other side of
                                ! the grid box, both in  [H L2 ~> m3 or kg].
  real, dimension(SZIB_(G)) :: &
    hlst, &             ! Work variable [H L2 ~> m3 or kg].
    Ihnew, &            ! Work variable [H-1 L-2 ~> m-3 or kg-1].
    CFL                 ! A nondimensional work variable.
  real :: min_h         ! The minimum thickness that can be realized during
                        ! any of the passes [H ~> m or kg m-2].
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
  logical :: do_j_tr(SZJ_(G))   ! If true, calculate the tracer profiles.
  logical :: do_i(SZIB_(G), SZJ_(G))     ! If true, work on given points.
  logical :: do_any_i
  integer :: i, j, j2, m, n, j_up, stencil
  real :: aR, aL, dMx, dMn, Tp, Tc, Tm, dA, mA, a6
  real :: fac1,v_L_in,v_L_out  ! terms used for time-stepping OBC reservoirs
  type(OBC_segment_type), pointer :: segment=>NULL()
  logical :: usePLMslope

  usePLMslope = .not. (usePPM .and. useHuynh)
  ! stencil for calculating slope values
  stencil = 1
  if (usePPM .and. .not. useHuynh) stencil = 2

  min_h = 0.1*GV%Angstrom_H
  h_neglect = GV%H_subroundoff
  !do i=is,ie ; ts2(i) = 0.0 ; enddo

  ! We conditionally perform work on tracer points: calculating the PLM slope,
  ! and updating tracer concentration within a cell
  ! this depends on whether there is a flux which would affect this tracer point,
  ! as indicated by domore_v. In the case of PPM reconstruction, a flux requires
  ! slope calculations at the two tracer points on either side (as indicated by
  ! the stencil variable), so we account for this with the do_j_tr flag array
  !
  ! Note: this does lead to unnecessary work in updating tracer concentrations,
  ! since that doesn't need a wider stencil with the PPM advection scheme, but
  ! this would require an additional loop, etc.
  do_j_tr(:) = .false.
  do J=js-1,je ; if (domore_v(J,k)) then ; do j2=1-stencil,stencil ; do_j_tr(j+j2) = .true. ; enddo ; endif ; enddo

  ! Calculate the j-direction profiles (slopes) of each tracer that
  ! is being advected.
  if (usePLMslope) then
    do j=js-stencil,je+stencil ; if (do_j_tr(j)) then ; do m=1,ntr ; do i=is,ie
      !if (ABS(Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j,k)) < &
      !    ABS(Tr(m)%t(i,j,k)-Tr(m)%t(i,j-1,k))) then
      !  maxslope = 4.0*(Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j,k))
      !else
      !  maxslope = 4.0*(Tr(m)%t(i,j,k)-Tr(m)%t(i,j-1,k))
      !endif
      !if ((Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j,k))*(Tr(m)%t(i,j,k)-Tr(m)%t(i,j-1,k)) < 0.0) then
      !  slope_y(i,m,j) = 0.0
      !elseif (ABS(Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j-1,k))<ABS(maxslope)) then
      !  slope_y(i,m,j) = G%mask2dCv(i,J) * G%mask2dCv(i,J-1) * &
      !                 0.5*(Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j-1,k))
      !else
      !  slope_y(i,m,j) = G%mask2dCv(i,J) * G%mask2dCv(i,J-1) * 0.5*maxslope
      !endif
       Tp = Tr(m)%t(i,j+1,k) ; Tc = Tr(m)%t(i,j,k) ; Tm = Tr(m)%t(i,j-1,k)
       dMx = max( Tp, Tc, Tm ) - Tc
       dMn= Tc - min( Tp, Tc, Tm )
       slope_y(i,m,j) = G%mask2dCv(i,J)*G%mask2dCv(i,J-1) * &
           sign( min(0.5*abs(Tp-Tm), 2.0*dMx, 2.0*dMn), Tp-Tm )
    enddo ; enddo ; endif ; enddo ! End of i-, m-, & j- loops.
  endif ! usePLMslope


  ! make a copy of the tracers in case values need to be overridden for OBCs

  do j=G%jsd,G%jed; do m=1,ntr; do i=G%isd,G%ied
    T_tmp(i,m,j) = Tr(m)%t(i,j,k)
  enddo ; enddo ; enddo

  ! loop through open boundaries and recalculate flux terms
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
     do n=1,OBC%number_of_segments
       segment=>OBC%segment(n)
       if (.not. associated(segment%tr_Reg)) cycle
       do i=is,ie
         if (segment%is_N_or_S) then
           if (i>=segment%HI%isd .and. i<=segment%HI%ied) then
              J = segment%HI%JsdB
              do m = 1,ntr ! replace tracers with OBC values
                if (associated(segment%tr_Reg%Tr(m)%tres)) then
                   if (segment%direction == OBC_DIRECTION_S) then
                      T_tmp(i,m,j) = segment%tr_Reg%Tr(m)%tres(i,j,k)
                   else
                      T_tmp(i,m,j+1) = segment%tr_Reg%Tr(m)%tres(i,j,k)
                   endif
                else
                   if (segment%direction == OBC_DIRECTION_S) then
                      T_tmp(i,m,j) = segment%tr_Reg%Tr(m)%OBC_inflow_conc
                   else
                      T_tmp(i,m,j+1) = segment%tr_Reg%Tr(m)%OBC_inflow_conc
                   endif
                endif
              enddo
              do m = 1,ntr ! Apply update tracer values for slope calculation
                do j=segment%HI%JsdB-1,segment%HI%JsdB+1
                  Tp = T_tmp(i,m,j+1) ; Tc = T_tmp(i,m,j) ; Tm = T_tmp(i,m,j-1)
                  dMx = max( Tp, Tc, Tm ) - Tc
                  dMn= Tc - min( Tp, Tc, Tm )
                  slope_y(i,m,j) = G%mask2dCv(i,J)*G%mask2dCv(i,J-1) * &
                       sign( min(0.5*abs(Tp-Tm), 2.0*dMx, 2.0*dMn), Tp-Tm )
                enddo
              enddo
           endif
         endif ! is_N_S
       enddo ! i-loop
     enddo ! segment loop
  endif; endif

  ! Calculate the j-direction fluxes of each tracer, using as much
  ! the minimum of the remaining mass flux (vhr) and the half the mass
  ! in the cell plus whatever part of its half of the mass flux that
  ! the flux through the other side does not require.
  do J=js-1,je ; if (domore_v(J,k)) then
    domore_v(J,k) = .false.

    do i=is,ie
      if (vhr(i,J,k) == 0.0) then
        vhh(i,J) = 0.0
        CFL(i) = 0.0
      elseif (vhr(i,J,k) < 0.0) then
        hup = hprev(i,j+1,k) - G%areaT(i,j+1)*min_h
        hlos = MAX(0.0,vhr(i,J+1,k))
        if ((((hup - hlos) + vhr(i,J,k)) < 0.0) .and. &
            ((0.5*hup + vhr(i,J,k)) < 0.0)) then
          vhh(i,J) = MIN(-0.5*hup,-hup+hlos,0.0)
          domore_v(J,k) = .true.
        else
          vhh(i,J) = vhr(i,J,k)
        endif
       !ts2(i) = 0.5*(1.0 + vhh(i,J) / (hprev(i,j+1,k) + h_neglect*G%areaT(i,j+1))
        CFL(i) = - vhh(i,J) / (hprev(i,j+1,k) + h_neglect*G%areaT(i,j+1)) ! CFL is positive
      else
        hup = hprev(i,j,k) - G%areaT(i,j)*min_h
        hlos = MAX(0.0,-vhr(i,J-1,k))
        if ((((hup - hlos) - vhr(i,J,k)) < 0.0) .and. &
            ((0.5*hup - vhr(i,J,k)) < 0.0)) then
          vhh(i,J) = MAX(0.5*hup,hup-hlos,0.0)
          domore_v(J,k) = .true.
        else
          vhh(i,J) = vhr(i,J,k)
        endif
       !ts2(i) = 0.5*(1.0 - vhh(i,J) / (hprev(i,j,k) + h_neglect*G%areaT(i,j)))
        CFL(i) = vhh(i,J) / (hprev(i,j,k) + h_neglect*G%areaT(i,j)) ! CFL is positive
      endif
    enddo

    if (usePPM) then
      do m=1,ntr ; do i=is,ie
        ! centre cell depending on upstream direction
        if (vhh(i,J) >= 0.0) then
          j_up = j
        else
          j_up = j + 1
        endif

        ! Implementation of PPM-H3
        Tp = Tr(m)%t(i,j_up+1,k) ; Tc = Tr(m)%t(i,j_up,k) ; Tm = Tr(m)%t(i,j_up-1,k)

        if (useHuynh) then
          aL = ( 5.*Tc + ( 2.*Tm - Tp ) )/6. ! H3 estimate
          aL = max( min(Tc,Tm), aL) ; aL = min( max(Tc,Tm), aL) ! Bound
          aR = ( 5.*Tc + ( 2.*Tp - Tm ) )/6. ! H3 estimate
          aR = max( min(Tc,Tp), aR) ; aR = min( max(Tc,Tp), aR) ! Bound
        else
          aL = 0.5 * ((Tm + Tc) + (slope_y(i,m,j_up-1) - slope_y(i,m,j_up)) / 3.)
          aR = 0.5 * ((Tc + Tp) + (slope_y(i,m,j_up) - slope_y(i,m,j_up+1)) / 3.)
        endif

        dA = aR - aL ; mA = 0.5*( aR + aL )
        if (G%mask2dCv(i,J_up)*G%mask2dCv(i,J_up-1)*(Tp-Tc)*(Tc-Tm) <= 0.) then
          aL = Tc ; aR = Tc ! PCM for local extremum and bounadry cells
        elseif ( dA*(Tc-mA) > (dA*dA)/6. ) then
          aL = 3.*Tc - 2.*aR
        elseif ( dA*(Tc-mA) < - (dA*dA)/6. ) then
          aR = 3.*Tc - 2.*aL
        endif

        a6 = 6.*Tc - 3. * (aR + aL) ! Curvature

        if (vhh(i,J) >= 0.0) then
          flux_y(i,m,J) = vhh(i,J)*( aR - 0.5 * CFL(i) * ( &
               ( aR - aL ) - a6 * ( 1. - 2./3. * CFL(I) ) ) )
        else
          flux_y(i,m,J) = vhh(i,J)*( aL + 0.5 * CFL(i) * ( &
               ( aR - aL ) + a6 * ( 1. - 2./3. * CFL(I) ) ) )
        endif
      enddo ; enddo
    else ! PLM
      do m=1,ntr ; do i=is,ie
        if (vhh(i,J) >= 0.0) then
          ! Indirect implementation of PLM
         !aL = Tr(m)%t(i,j,k) - 0.5 * slope_y(i,m,j)
         !aR = Tr(m)%t(i,j,k) + 0.5 * slope_y(i,m,j)
         !flux_y(i,m,J) = vhh(i,J)*( aR - 0.5 * (aR-aL) * CFL(i) )
          ! Alternative implementation of PLM
         !aR = Tr(m)%t(i,j,k) + 0.5 * slope_y(i,m,j)
         !flux_y(i,m,J) = vhh(i,J)*(aR - 0.5 * slope_y(i,m,j)*CFL(i))
          ! Alternative implementation of PLM
          Tc = Tr(m)%t(i,j,k)
          flux_y(i,m,J) = vhh(i,J)*( Tc + 0.5 * slope_y(i,m,j) * ( 1. - CFL(i) ) )
          ! Original implementation of PLM
         !flux_y(i,m,J) = vhh(i,J)*(Tr(m)%t(i,j,k) + slope_y(i,m,j)*ts2(i))
        else
          ! Indirect implementation of PLM
         !aL = Tr(m)%t(i,j+1,k) - 0.5 * slope_y(i,m,j+1)
         !aR = Tr(m)%t(i,j+1,k) + 0.5 * slope_y(i,m,j+1)
         !flux_y(i,m,J) = vhh(i,J)*( aL + 0.5 * (aR-aL) * CFL(i) )
          ! Alternative implementation of PLM
         !aL = Tr(m)%t(i,j+1,k) - 0.5 * slope_y(i,m,j+1)
         !flux_y(i,m,J) = vhh(i,J)*( aL + 0.5 * slope_y(i,m,j+1)*CFL(i) )
          ! Alternative implementation of PLM
          Tc = Tr(m)%t(i,j+1,k)
          flux_y(i,m,J) = vhh(i,J)*( Tc - 0.5 * slope_y(i,m,j+1) * ( 1. - CFL(i) ) )
          ! Original implementation of PLM
         !flux_y(i,m,J) = vhh(i,J)*(Tr(m)%t(i,j+1,k) - slope_y(i,m,j+1)*ts2(i))
        endif
      enddo ; enddo
    endif ! usePPM

    if (associated(OBC)) then ; if (OBC%OBC_pe) then
      if (OBC%specified_v_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally) then
        do n=1,OBC%number_of_segments
          segment=>OBC%segment(n)
          if (.not. segment%specified) cycle
          if (.not. associated(segment%tr_Reg)) cycle
          if (OBC%segment(n)%is_N_or_S) then
            if (J >= segment%HI%JsdB .and. J<= segment%HI%JedB) then
              do i=segment%HI%isd,segment%HI%ied
                ! Tracer fluxes are set to prescribed values only for inflows from masked areas.
                ! Now changing to simply fixed inflows.
                if ((vhr(i,J,k) > 0.0) .and. (segment%direction == OBC_DIRECTION_S) .or. &
                   (vhr(i,J,k) < 0.0) .and. (segment%direction == OBC_DIRECTION_N)) then
                  vhh(i,J) = vhr(i,J,k)
                  do m=1,ntr
                    if (associated(segment%tr_Reg%Tr(m)%t)) then
                      flux_y(i,m,J) = vhh(i,J)*OBC%segment(n)%tr_Reg%Tr(m)%tres(i,J,k)
                    else ; flux_y(i,m,J) = vhh(i,J)*OBC%segment(n)%tr_Reg%Tr(m)%OBC_inflow_conc ; endif
                  enddo
                endif
              enddo
            endif
          endif
        enddo
      endif

      if (OBC%open_v_BCs_exist_globally) then
        do n=1,OBC%number_of_segments
          segment=>OBC%segment(n)
          if (segment%specified) cycle
          if (.not. associated(segment%tr_Reg)) cycle
          if (segment%is_N_or_S .and. (J >= segment%HI%JsdB .and. J<= segment%HI%JedB)) then
            do i=segment%HI%isd,segment%HI%ied
              ! Tracer fluxes are set to prescribed values only for inflows from masked areas.
              if ((vhr(i,J,k) > 0.0) .and. (G%mask2dT(i,j) < 0.5) .or. &
                  (vhr(i,J,k) < 0.0) .and. (G%mask2dT(i,j+1) < 0.5)) then
                vhh(i,J) = vhr(i,J,k)
                do m=1,ntr
                  if (associated(segment%tr_Reg%Tr(m)%t)) then
                    flux_y(i,m,J) = vhh(i,J)*segment%tr_Reg%Tr(m)%tres(i,J,k)
                  else ; flux_y(i,m,J) = vhh(i,J)*segment%tr_Reg%Tr(m)%OBC_inflow_conc ; endif
                enddo
              endif
            enddo
          endif
        enddo
      endif
    endif; endif

  else ! not domore_v.
    do i=is,ie ; vhh(i,J) = 0.0 ; enddo
    do m=1,ntr ; do i=is,ie ; flux_y(i,m,J) = 0.0 ; enddo ; enddo
  endif ; enddo ! End of j-loop

  do J=js-1,je ; do i=is,ie
    vhr(i,J,k) = vhr(i,J,k) - vhh(i,J)
    if (abs(vhr(i,J,k)) < vh_neglect(i,J)) vhr(i,J,k) = 0.0
  enddo ; enddo

  ! Calculate new tracer concentration in each cell after accounting
  ! for the j-direction fluxes.
  do j=js,je ; if (do_j_tr(j)) then
    do i=is,ie
      if ((vhh(i,J) /= 0.0) .or. (vhh(i,J-1) /= 0.0)) then
        do_i(i,j) = .true.
        hlst(i) = hprev(i,j,k)
        hprev(i,j,k) = max(hprev(i,j,k) - (vhh(i,J) - vhh(i,J-1)), 0.0)
        if (hprev(i,j,k) <= 0.0) then ; do_i(i,j) = .false.
        elseif (hprev(i,j,k) < h_neglect*G%areaT(i,j)) then
          hlst(i) = hlst(i) + (h_neglect*G%areaT(i,j) - hprev(i,j,k))
          Ihnew(i) = 1.0 / (h_neglect*G%areaT(i,j))
        else ;  Ihnew(i) = 1.0 / hprev(i,j,k) ; endif
      else ; do_i(i,j) = .false. ; endif
    enddo

    ! update tracer and save some diagnostics
    do m=1,ntr
      do i=is,ie ; if (do_i(i,j)) then
        Tr(m)%t(i,j,k) = (Tr(m)%t(i,j,k) * hlst(i) - &
                          (flux_y(i,m,J) - flux_y(i,m,J-1))) * Ihnew(i)
      endif ; enddo

      ! diagnostics
      if (associated(Tr(m)%ad_y)) then ; do i=is,ie ; if (do_i(i,j)) then
        Tr(m)%ad_y(i,J,k) = Tr(m)%ad_y(i,J,k) + flux_y(i,m,J)*Idt
      endif ; enddo ; endif

      ! diagnose convergence of flux_y and add to convergence of flux_x.
      ! division by areaT to get into W/m2 for heat and kg/(s*m2) for salt.
      if (associated(Tr(m)%advection_xy)) then
        do i=is,ie ; if (do_i(i,j)) then
          Tr(m)%advection_xy(i,j,k) = Tr(m)%advection_xy(i,j,k) - (flux_y(i,m,J) - flux_y(i,m,J-1))* Idt * &
                                          G%IareaT(i,j)
        endif ; enddo
      endif

    enddo
  endif ; enddo ! End of j-loop.

  ! compute ad2d_y diagnostic outside above j-loop so as to make the summation ordered when OMP is active.

  !$OMP ordered
  do j=js,je ; if (do_j_tr(j)) then
    do m=1,ntr
      if (associated(Tr(m)%ad2d_y)) then ; do i=is,ie ; if (do_i(i,j)) then
        Tr(m)%ad2d_y(i,J) = Tr(m)%ad2d_y(i,J) + flux_y(i,m,J)*Idt
      endif ; enddo ; endif
    enddo
  endif ; enddo ! End of j-loop.
  !$OMP end ordered

end subroutine advect_y

!> Initialize lateral tracer advection module
subroutine tracer_advect_init(Time, G, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time        !< current model time
  type(ocean_grid_type),   intent(in)    :: G           !< ocean grid structure
  type(unit_scale_type),   intent(in)    :: US          !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file  !< open file to parse for model parameters
  type(diag_ctrl), target, intent(inout) :: diag        !< regulates diagnostic output
  type(tracer_advect_CS),  pointer       :: CS          !< module control structure

  integer, save :: init_calls = 0

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_tracer_advect" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (associated(CS)) then
    call MOM_error(WARNING, "tracer_advect_init called with associated control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DT", CS%dt, fail_if_missing=.true., &
          desc="The (baroclinic) dynamics time step.", units="s", scale=US%s_to_T)
  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false.)
  call get_param(param_file, mdl, "TRACER_ADVECTION_SCHEME", mesg, &
          desc="The horizontal transport scheme for tracers:\n"//&
          "  PLM    - Piecewise Linear Method\n"//&
          "  PPM:H3 - Piecewise Parabolic Method (Huyhn 3rd order)\n"// &
          "  PPM    - Piecewise Parabolic Method (Colella-Woodward)" &
          , default='PLM')
  select case (trim(mesg))
    case ("PLM")
      CS%usePPM = .false.
    case ("PPM:H3")
      CS%usePPM = .true.
      CS%useHuynh = .true.
    case ("PPM")
      CS%usePPM = .true.
      CS%useHuynh = .false.
    case default
      call MOM_error(FATAL, "MOM_tracer_advect, tracer_advect_init: "//&
           "Unknown TRACER_ADVECTION_SCHEME = "//trim(mesg))
  end select

  id_clock_advect = cpu_clock_id('(Ocean advect tracer)', grain=CLOCK_MODULE)
  id_clock_pass = cpu_clock_id('(Ocean tracer halo updates)', grain=CLOCK_ROUTINE)
  id_clock_sync = cpu_clock_id('(Ocean tracer global synch)', grain=CLOCK_ROUTINE)

end subroutine tracer_advect_init

!> Close the tracer advection module
subroutine tracer_advect_end(CS)
  type(tracer_advect_CS), pointer :: CS  !< module control structure

  if (associated(CS)) deallocate(CS)

end subroutine tracer_advect_end


!> \namespace mom_tracer_advect
!!
!!    This program contains the subroutines that advect tracers
!!  horizontally (i.e. along layers).
!!
!! \section section_mom_advect_intro
!!
!!  * advect_tracer advects tracer concentrations using a combination
!!  of the modified flux advection scheme from Easter (Mon. Wea. Rev.,
!!  1993) with tracer distributions given by the monotonic modified
!!  van Leer scheme proposed by Lin et al. (Mon. Wea. Rev., 1994).
!!  This scheme conserves the total amount of tracer while avoiding
!!  spurious maxima and minima of the tracer concentration.  If a
!!  higher order accuracy scheme is needed, suggest monotonic
!!  piecewise parabolic method, as described in Carpenter et al.
!!  (MWR, 1990).
!!
!!  * advect_tracer has 4 arguments, described below. This
!!  subroutine determines the volume of a layer in a grid cell at the
!!  previous instance when the tracer concentration was changed, so
!!  it is essential that the volume fluxes should be correct.  It is
!!  also important that the tracer advection occurs before each
!!  calculation of the diabatic forcing.

end module MOM_tracer_advect
