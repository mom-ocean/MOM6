!>  This program contains the subroutines that advect tracers
!!  along coordinate surfaces.
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
use MOM_tracer_registry, only : tracer_registry_type, tracer_type
use MOM_verticalGrid,    only : verticalGrid_type
implicit none ; private

#include <MOM_memory.h>

public advect_tracer
public tracer_advect_init
public tracer_advect_end

!> Control structure for this module
type, public :: tracer_advect_CS ; private
  real    :: dt                    !< The baroclinic dynamics time step, in s.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !< timing of diagnostic output.
  logical :: debug                 !< If true, write verbose checksums for debugging purposes.
  logical :: usePPM                !< If true, use PPM instead of PLM
  logical :: useHuynh              !< If true, use the Huynh scheme for PPM interface values
  type(group_pass_type) :: pass_uhr_vhr_t_hprev ! For group pass
end type tracer_advect_CS

integer :: id_clock_advect
integer :: id_clock_pass
integer :: id_clock_sync

contains

!> This routine time steps the tracer concentration using a
!! monotonic, conservative, weakly diffusive scheme.
subroutine advect_tracer(h_end, uhtr, vhtr, OBC, dt, G, GV, CS, Reg, &
      h_prev_opt, max_iter_in, x_first_in, uhr_out, vhr_out, h_out)
  type(ocean_grid_type),                     intent(inout) :: G     !< ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV    !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h_end !< layer thickness after advection (m or kg m-2)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: uhtr  !< accumulated volume/mass flux through zonal face (m3 or kg)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: vhtr  !< accumulated volume/mass flux through merid face (m3 or kg)
  type(ocean_OBC_type),                      pointer       :: OBC   !< specifies whether, where, and what OBCs are used
  real,                                      intent(in)    :: dt    !< time increment (seconds)
  type(tracer_advect_CS),                    pointer       :: CS    !< control structure for module
  type(tracer_registry_type),                pointer       :: Reg   !< pointer to tracer registry
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  optional      :: h_prev_opt !< layer thickness before advection (m or kg m-2)
  integer,                                   optional      :: max_iter_in
  logical,                                   optional      :: x_first_in
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), optional, intent(out)    :: uhr_out  !< accumulated volume/mass flux through zonal face (m3 or kg)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), optional, intent(out)    :: vhr_out  !< accumulated volume/mass flux through merid face (m3 or kg)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  optional      :: h_out !< layer thickness before advection (m or kg m-2)

  type(tracer_type) :: Tr(MAX_FIELDS_) ! The array of registered tracers
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    hprev           ! cell volume at the end of previous tracer change (m3)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: &
    uhr             ! The remaining zonal thickness flux (m3)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: &
    vhr             ! The remaining meridional thickness fluxes (m3)
  real :: uh_neglect(SZIB_(G),SZJ_(G)) ! uh_neglect and vh_neglect are the
  real :: vh_neglect(SZI_(G),SZJB_(G)) ! magnitude of remaining transports that
                                       ! can be simply discarded (m3 or kg).

  real :: landvolfill                   ! An arbitrary? nonzero cell volume, m3.
  real :: Idt                           ! 1/dt in s-1.
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

  ! increase stencil size for Colella & Woodward PPM
  if (CS%usePPM .and. .not. CS%useHuynh) stencil = 3

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_tracer_advect: "// &
       "tracer_advect_init must be called before advect_tracer.")
  if (.not. associated(Reg)) call MOM_error(FATAL, "MOM_tracer_advect: "// &
       "register_tracer must be called before advect_tracer.")
  if (Reg%ntr==0) return
  call cpu_clock_begin(id_clock_advect)
  x_first = (MOD(G%first_direction,2) == 0)

  ntr = Reg%ntr
  do m=1,ntr ; Tr(m) = Reg%Tr(m) ; enddo
  Idt = 1.0/dt

  max_iter = 2*INT(CEILING(dt/CS%dt)) + 1

  if(present(max_iter_in)) max_iter = max_iter_in
  if(present(x_first_in))  x_first = x_first_in
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

  do k = 1, nz
    do j = jsd,  jed;  do i = IsdB, IedB; uhr(i,j,k) = 0.0; enddo ; enddo
    do j = jsdB, jedB; do i = Isd,  Ied;  vhr(i,j,k) = 0.0; enddo ; enddo
    do j = jsd,  jed;  do i = Isd,  Ied;  hprev(i,j,k) = 0.0; enddo ; enddo
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
        hprev(i,j,k) = h_prev_opt(i,j,k);
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
      call cpu_clock_begin(id_clock_pass)
      call do_group_pass(CS%pass_uhr_vhr_t_hprev, G%Domain)
      call cpu_clock_end(id_clock_pass)

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

!$OMP parallel do default(none) shared(nz,domore_k,x_first,Tr,hprev,uhr,uh_neglect,  &
!$OMP                                  OBC,domore_u,ntr,Idt,isv,iev,jsv,jev,stencil, &
!$OMP                                  G,GV,CS,vhr,vh_neglect,domore_v)

    !  To ensure positive definiteness of the thickness at each iteration, the
    !  mass fluxes out of each layer are checked each step, and limited to keep
    !  the thicknesses positive.  This means that several iterations may be required
    !  for all the transport to happen.  The sum over domore_k keeps the processors
    !  synchronized.  This may not be very efficient, but it should be reliable.
    do k=1,nz ; if (domore_k(k) > 0) then

      if (x_first) then

        ! First, advect zonally.
        call advect_x(Tr, hprev, uhr, uh_neglect, OBC, domore_u, ntr, Idt, &
                      isv, iev, jsv-stencil, jev+stencil, k, G, GV, CS%usePPM, CS%useHuynh)

        !  Next, advect meridionally.
        call advect_y(Tr, hprev, vhr, vh_neglect, OBC, domore_v, ntr, Idt, &
                      isv, iev, jsv, jev, k, G, GV, CS%usePPM, CS%useHuynh)

        domore_k(k) = 0
        do j=jsv-stencil,jev+stencil ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
        do J=jsv-1,jev ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo

      else

        ! First, advect meridionally.
        call advect_y(Tr, hprev, vhr, vh_neglect, OBC, domore_v, ntr, Idt, &
                      isv-stencil, iev+stencil, jsv, jev, k, G, GV, CS%usePPM, CS%useHuynh)

        ! Next, advect zonally.
        call advect_x(Tr, hprev, uhr, uh_neglect, OBC, domore_u, ntr, Idt, &
                      isv, iev, jsv, jev, k, G, GV, CS%usePPM, CS%useHuynh)

        domore_k(k) = 0
        do j=jsv,jev ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
        do J=jsv-1,jev ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo

      endif


    endif ; enddo ! End of k-loop

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

  if(present(uhr_out)) uhr_out(:,:,:) = uhr(:,:,:)
  if(present(vhr_out)) vhr_out(:,:,:) = vhr(:,:,:)
  if(present(h_out)) h_out(:,:,:) = hprev(:,:,:)

  call cpu_clock_end(id_clock_advect)

end subroutine advect_tracer


!> This subroutine does 1-d flux-form advection in the zonal direction using
!! a monotonic piecewise linear scheme.
subroutine advect_x(Tr, hprev, uhr, uh_neglect, OBC, domore_u, ntr, Idt, &
                    is, ie, js, je, k, G, GV, usePPM, useHuynh)
  type(ocean_grid_type),                     intent(inout) :: G
  type(verticalGrid_type),                   intent(in)    :: GV
  type(tracer_type), dimension(ntr),         intent(inout) :: Tr
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: hprev
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uhr
   real, dimension(SZIB_(G),SZJ_(G)),        intent(inout) :: uh_neglect
  type(ocean_OBC_type),                      pointer       :: OBC
  logical, dimension(SZJ_(G),SZK_(G)),       intent(inout) :: domore_u
  real,                                      intent(in)    :: Idt
  integer,                                   intent(in)    :: ntr, is, ie, js, je,k
  logical,                                   intent(in)    :: usePPM, useHuynh

  real, dimension(SZIB_(G),ntr) :: &
    slope_x, &          ! The concentration slope per grid point in units of
                        ! concentration (nondim.).
    flux_x              ! The tracer flux across a boundary in m3*conc or kg*conc.
  real :: maxslope      ! The maximum concentration slope per grid point
                        ! consistent with monotonicity, in conc. (nondim.).
  real :: hup, hlos     ! hup is the upwind volume, hlos is the
                        ! part of that volume that might be lost
                        ! due to advection out the other side of
                        ! the grid box, both in m3 or kg.
  real :: uhh(SZIB_(G)) ! The zonal flux that occurs during the
                        ! current iteration, in m3 or kg.
  real, dimension(SZIB_(G)) :: &
    hlst, Ihnew, &      ! Work variables with units of m3 or kg and m-3 or kg-1.
    CFL                 ! A nondimensional work variable.
  real :: min_h         ! The minimum thickness that can be realized during
                        ! any of the passes, in m or kg m-2.
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected, in m.
  logical :: do_i(SZIB_(G))     ! If true, work on given points.
  logical :: do_any_i
  integer :: i, j, m, n, i_up
  real :: aR, aL, dMx, dMn, Tp, Tc, Tm, dA, mA, a6
  logical :: usePLMslope

  usePLMslope = .not. (usePPM .and. useHuynh)

  min_h = 0.1*GV%Angstrom
  h_neglect = GV%H_subroundoff

! do I=is-1,ie ; ts2(I) = 0.0 ; enddo
  do I=is-1,ie ; CFL(I) = 0.0 ; enddo

  do j=js,je ; if (domore_u(j,k)) then
    domore_u(j,k) = .false.

    ! Calculate the i-direction profiles (slopes) of each tracer that
    ! is being advected.
    if (usePLMslope) then
      do m=1,ntr ; do i=is-1,ie+1
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
       !ts2(I) = 0.5*(1.0 + uhh(I)/(hprev(i+1,j,k)+h_neglect))
        CFL(I) = - uhh(I)/(hprev(i+1,j,k)+h_neglect) ! CFL is positive
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
       !ts2(I) = 0.5*(1.0 - uhh(I)/(hprev(i,j,k)+h_neglect))
        CFL(I) = uhh(I)/(hprev(i,j,k)+h_neglect) ! CFL is positive
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
        Tp = Tr(m)%t(i_up+1,j,k) ; Tc = Tr(m)%t(i_up,j,k) ; Tm = Tr(m)%t(i_up-1,j,k)

        if (useHuynh) then
          aL = ( 5.*Tc + ( 2.*Tm - Tp ) )/6. ! H3 estimate
          aL = max( min(Tc,Tm), aL) ; aL = min( max(Tc,Tm), aL) ! Bound
          aR = ( 5.*Tc + ( 2.*Tp - Tm ) )/6. ! H3 estimate
          aR = max( min(Tc,Tp), aR) ; aR = min( max(Tc,Tp), aR) ! Bound
        else
          aL = 0.5 * ((Tm + Tc) + (slope_x(i_up-1,m) - slope_x(i_up,m)))
          aR = 0.5 * ((Tc + Tp) + (slope_x(i_up,m) - slope_x(i_up+1,m)))
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
          flux_x(I,m) = uhh(I)*( aR - 0.5 * CFL(I) * ( &
               ( aR - aL ) - a6 * ( 1. - 2./3. * CFL(I) ) ) )
        else
          flux_x(I,m) = uhh(I)*( aL + 0.5 * CFL(I) * ( &
               ( aR - aL ) + a6 * ( 1. - 2./3. * CFL(I) ) ) )
        endif
      enddo ; enddo
    else ! PLM
      do m=1,ntr ; do I=is-1,ie
        if (uhh(I) >= 0.0) then
          ! Indirect implementation of PLM
         !aL = Tr(m)%t(i,j,k) - 0.5 * slope_x(i,m)
         !aR = Tr(m)%t(i,j,k) + 0.5 * slope_x(i,m)
         !flux_x(I,m) = uhh(I)*( aR - 0.5 * (aR-aL) * CFL(I) )
          ! Alternative implementation of PLM
         !aR = Tr(m)%t(i,j,k) + 0.5 * slope_x(i,m)
         !flux_x(I,m) = uhh(I)*( aR - 0.5 * slope_x(i,m) * CFL(I) )
          ! Alternative implementation of PLM
          Tc = Tr(m)%t(i,j,k)
          flux_x(I,m) = uhh(I)*( Tc + 0.5 * slope_x(i,m) * ( 1. - CFL(I) ) )
          ! Original implementation of PLM
         !flux_x(I,m) = uhh(I)*(Tr(m)%t(i,j,k) + slope_x(i,m)*ts2(I))
        else
          ! Indirect implementation of PLM
         !aL = Tr(m)%t(i+1,j,k) - 0.5 * slope_x(i+1,m)
         !aR = Tr(m)%t(i+1,j,k) + 0.5 * slope_x(i+1,m)
         !flux_x(I,m) = uhh(I)*( aL + 0.5 * (aR-aL) * CFL(I) )
          ! Alternative implementation of PLM
         !aL = Tr(m)%t(i+1,j,k) - 0.5 * slope_x(i+1,m)
         !flux_x(I,m) = uhh(I)*( aL + 0.5 * slope_x(i+1,m) * CFL(I) )
          ! Alternative implementation of PLM
          Tc = Tr(m)%t(i+1,j,k)
          flux_x(I,m) = uhh(I)*( Tc - 0.5 * slope_x(i+1,m) * ( 1. - CFL(I) ) )
          ! Original implementation of PLM
         !flux_x(I,m) = uhh(I)*(Tr(m)%t(i+1,j,k) - slope_x(i+1,m)*ts2(I))
        endif
       !ts2(I) = 0.5*(1.0 - uhh(I)/(hprev(i,j,k)+h_neglect))
      enddo ; enddo
    endif ! usePPM

    if (associated(OBC)) then ; if (OBC%OBC_pe) then ; if (OBC%specified_u_BCs_exist_globally) then
      do n=1,OBC%number_of_segments
        if (OBC%segment(n)%is_E_or_W) then
          if (j >= OBC%segment(n)%HI%jsd .and. j<= OBC%segment(n)%HI%jed) then
            I = OBC%segment(n)%HI%IsdB
            ! Tracer fluxes are set to prescribed values only for inflows from masked areas.
            if ((uhr(I,j,k) > 0.0) .and. (G%mask2dT(i,j) < 0.5) .or. &
                (uhr(I,j,k) < 0.0) .and. (G%mask2dT(i+1,j) < 0.5)) then
              uhh(I) = uhr(I,j,k)
              do m=1,ntr
                if (associated(Tr(m)%OBC_in_u)) then
                  flux_x(I,m) = uhh(I)*Tr(m)%OBC_in_u(I,j,k)
                else ; flux_x(I,m) = uhh(I)*Tr(m)%OBC_inflow_conc ; endif
              enddo
            endif
          endif
        endif
      enddo
    endif ; endif ; endif

    ! Calculate new tracer concentration in each cell after accounting
    ! for the i-direction fluxes.
    do I=is-1,ie
      uhr(I,j,k) = uhr(I,j,k) - uhh(I)
      if (abs(uhr(I,j,k)) < uh_neglect(I,j)) uhr(I,j,k) = 0.0
    enddo
    do i=is,ie
      if ((uhh(I) /= 0.0) .or. (uhh(I-1) /= 0.0)) then
        do_i(i) = .true.
        hlst(i) = hprev(i,j,k)
        hprev(i,j,k) = hprev(i,j,k) - (uhh(I) - uhh(I-1))
        if (hprev(i,j,k) <= 0.0) then ; do_i(i) = .false.
        elseif (hprev(i,j,k) < h_neglect*G%areaT(i,j)) then
          hlst(i) = hlst(i) + (h_neglect*G%areaT(i,j) - hprev(i,j,k))
          Ihnew(i) = 1.0 / (h_neglect*G%areaT(i,j))
        else ;  Ihnew(i) = 1.0 / hprev(i,j,k) ; endif
      else
        do_i(i) = .false.
      endif
    enddo

    ! update tracer concentration from i-flux and save some diagnostics
    do m=1,ntr

      ! update tracer
      do i=is,ie ; if ((do_i(i)) .and. (Ihnew(i) > 0.0)) then
        Tr(m)%t(i,j,k) = (Tr(m)%t(i,j,k) * hlst(i) - &
                          (flux_x(I,m) - flux_x(I-1,m))) * Ihnew(i)
      endif ; enddo

      ! diagnostics
      if (associated(Tr(m)%ad_x)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad_x(I,j,k) = Tr(m)%ad_x(I,j,k) + flux_x(I,m)*Idt
      endif ; enddo ; endif
      if (associated(Tr(m)%ad2d_x)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad2d_x(I,j) = Tr(m)%ad2d_x(I,j) + flux_x(I,m)*Idt
      endif ; enddo ; endif

      ! diagnose convergence of flux_x (do not use the Ihnew(i) part of the logic).
      ! division by areaT to get into W/m2 for heat and kg/(s*m2) for salt.
      if (associated(Tr(m)%advection_xy)) then
        do i=is,ie ; if (do_i(i)) then
          Tr(m)%advection_xy(i,j,k) = Tr(m)%advection_xy(i,j,k) - (flux_x(I,m) - flux_x(I-1,m)) * Idt * G%IareaT(i,j)
        endif ; enddo
      endif

    enddo

  endif ; enddo ! End of j-loop.

end subroutine advect_x

!> This subroutine does 1-d flux-form advection using a monotonic piecewise
!! linear scheme.
subroutine advect_y(Tr, hprev, vhr, vh_neglect, OBC, domore_v, ntr, Idt, &
                    is, ie, js, je, k, G, GV, usePPM, useHuynh)
  type(ocean_grid_type),                     intent(inout) :: G
  type(verticalGrid_type),                   intent(in)    :: GV
  type(tracer_type), dimension(ntr),         intent(inout) :: Tr
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: hprev
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vhr
  real, dimension(SZI_(G),SZJB_(G)),         intent(inout) :: vh_neglect
  type(ocean_OBC_type),                      pointer       :: OBC
  logical, dimension(SZJB_(G),SZK_(G)),      intent(inout) :: domore_v
  real,                                      intent(in)    :: Idt
  integer,                                   intent(in)    :: ntr, is, ie, js, je,k
  logical,                                   intent(in)    :: usePPM, useHuynh

  real, dimension(SZI_(G),ntr,SZJB_(G)) :: &
    slope_y, &                  ! The concentration slope per grid point in units of
                                ! concentration (nondim.).
    flux_y                      ! The tracer flux across a boundary in m3 * conc or kg*conc.
  real :: maxslope              ! The maximum concentration slope per grid point
                                ! consistent with monotonicity, in conc. (nondim.).
  real :: vhh(SZI_(G),SZJB_(G)) ! The meridional flux that occurs during the
                                ! current iteration, in m3 or kg.
  real :: hup, hlos             ! hup is the upwind volume, hlos is the
                                ! part of that volume that might be lost
                                ! due to advection out the other side of
                                ! the grid box, both in m3 or kg.
  real, dimension(SZIB_(G)) :: &
    hlst, Ihnew, &      ! Work variables with units of m3 or kg and m-3 or kg-1.
    CFL                 ! A nondimensional work variable.
  real :: min_h         ! The minimum thickness that can be realized during
                        ! any of the passes, in m or kg m-2.
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected, in m.
  logical :: do_j_tr(SZJ_(G))   ! If true, calculate the tracer profiles.
  logical :: do_i(SZIB_(G))     ! If true, work on given points.
  logical :: do_any_i
  integer :: i, j, m, n, j_up
  real :: aR, aL, dMx, dMn, Tp, Tc, Tm, dA, mA, a6
  logical :: usePLMslope

  usePLMslope = .not. (usePPM .and. useHuynh)

  min_h = 0.1*GV%Angstrom
  h_neglect = GV%H_subroundoff

 !do i=is,ie ; ts2(i) = 0.0 ; enddo
  do_j_tr(js-1) = domore_v(js-1,k) ; do_j_tr(je+1) = domore_v(je,k)
  do j=js,je ; do_j_tr(j) = (domore_v(J-1,k) .or. domore_v(J,k)) ; enddo

  !   Calculate the j-direction profiles (slopes) of each tracer that
  ! is being advected.
  if (usePLMslope) then
    do j=js-1,je+1 ; if (do_j_tr(j)) then ; do m=1,ntr ; do i=is,ie
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
       !ts2(i) = 0.5*(1.0 + vhh(i,J) / (hprev(i,j+1,k)+h_neglect))
        CFL(i) = - vhh(i,J) / (hprev(i,j+1,k)+h_neglect) ! CFL is positive
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
       !ts2(i) = 0.5*(1.0 - vhh(i,J) / (hprev(i,j,k)+h_neglect))
        CFL(i) = vhh(i,J) / (hprev(i,j,k)+h_neglect) ! CFL is positive
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
          aL = 0.5 * ((Tm + Tc) + (slope_y(i,m,j_up-1) - slope_y(i,m,j_up)))
          aR = 0.5 * ((Tc + Tp) + (slope_y(i,m,j_up) - slope_y(i,m,j_up+1)))
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

    if (associated(OBC)) then ; if (OBC%OBC_pe) then ; if (OBC%specified_v_BCs_exist_globally) then
      do n=1,OBC%number_of_segments
        if (OBC%segment(n)%is_N_or_S) then
          if (J >= OBC%segment(n)%HI%JsdB .and. J<= OBC%segment(n)%HI%JedB) then
            i = OBC%segment(n)%HI%isd
            ! Tracer fluxes are set to prescribed values only for inflows from masked areas.
            if ((vhr(i,J,k) > 0.0) .and. (G%mask2dT(i,j) < 0.5) .or. &
                (vhr(i,J,k) < 0.0) .and. (G%mask2dT(i,j+1) < 0.5)) then
              vhh(i,J) = vhr(i,J,k)
              do m=1,ntr
                if (associated(Tr(m)%OBC_in_v)) then
                  flux_y(i,m,J) = vhh(i,J)*Tr(m)%OBC_in_v(i,J,k)
                else ; flux_y(i,m,J) = vhh(i,J)*Tr(m)%OBC_inflow_conc ; endif
              enddo
            endif
          endif
        endif
      enddo
    endif ; endif ; endif
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
        do_i(i) = .true.
        hlst(i) = hprev(i,j,k)
        hprev(i,j,k) = max(hprev(i,j,k) - (vhh(i,J) - vhh(i,J-1)), 0.0)
        if (hprev(i,j,k) <= 0.0) then ; do_i(i) = .false.
        elseif (hprev(i,j,k) < h_neglect*G%areaT(i,j)) then
          hlst(i) = hlst(i) + (h_neglect*G%areaT(i,j) - hprev(i,j,k))
          Ihnew(i) = 1.0 / (h_neglect*G%areaT(i,j))
        else ;  Ihnew(i) = 1.0 / hprev(i,j,k) ; endif
      else ; do_i(i) = .false. ; endif
    enddo

    ! update tracer and save some diagnostics
    do m=1,ntr
      do i=is,ie ; if (do_i(i)) then
        Tr(m)%t(i,j,k) = (Tr(m)%t(i,j,k) * hlst(i) - &
                          (flux_y(i,m,J) - flux_y(i,m,J-1))) * Ihnew(i)
      endif ; enddo

      ! diagnostics
      if (associated(Tr(m)%ad_y)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad_y(i,J,k) = Tr(m)%ad_y(i,J,k) + flux_y(i,m,J)*Idt
      endif ; enddo ; endif
      if (associated(Tr(m)%ad2d_y)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad2d_y(i,J) = Tr(m)%ad2d_y(i,J) + flux_y(i,m,J)*Idt
      endif ; enddo ; endif

      ! diagnose convergence of flux_y and add to convergence of flux_x.
      ! division by areaT to get into W/m2 for heat and kg/(s*m2) for salt.
      if (associated(Tr(m)%advection_xy)) then
        do i=is,ie ; if (do_i(i)) then
          Tr(m)%advection_xy(i,j,k) = Tr(m)%advection_xy(i,j,k) - (flux_y(i,m,J) - flux_y(i,m,J-1))* Idt * G%IareaT(i,j)
        endif ; enddo
      endif


    enddo
  endif ; enddo ! End of j-loop.


end subroutine advect_y

!> Initialize lateral tracer advection module
subroutine tracer_advect_init(Time, G, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time        !< current model time
  type(ocean_grid_type),   intent(in)    :: G           !< ocean grid structure
  type(param_file_type),   intent(in)    :: param_file  !< open file to parse for model parameters
  type(diag_ctrl), target, intent(inout) :: diag        !< regulates diagnostic output
  type(tracer_advect_CS),  pointer       :: CS          !< module control structure

  integer, save :: init_calls = 0

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_tracer_advect" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (associated(CS)) then
    call MOM_error(WARNING, "tracer_advect_init called with associated control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "DT", CS%dt, fail_if_missing=.true., &
          desc="The (baroclinic) dynamics time step.", units="s")
  call get_param(param_file, mod, "DEBUG", CS%debug, default=.false.)
  call get_param(param_file, mod, "TRACER_ADVECTION_SCHEME", mesg, &
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
  type(tracer_advect_CS), pointer :: CS

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
!!
!!  \section section_gridlayout MOM grid layout
!!
!!  A small fragment of the grid is shown below:
!!
!! \verbatim
!!    j+1  x ^ x ^ x
!!
!!    j+1  > o > o >
!!
!!    j    x ^ x ^ x
!!
!!    j    > o > o >
!!
!!    j-1  x ^ x ^ x
!!
!!        i-1  i  i+1
!!
!!           i  i+1
!!
!! \endverbatim
!!
!!  Fields at each point
!!  * x =  q, CoriolisBu
!!  * ^ =  v, PFv, CAv, vh, diffv, tauy, vbt, vhtr
!!  * > =  u, PFu, CAu, uh, diffu, taux, ubt, uhtr
!!  * o =  h, bathyT, eta, T, S, tr
!!
!!  The boundaries always run through q grid points (x).
!!

end module MOM_tracer_advect
