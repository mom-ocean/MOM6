module MOM_tracer
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, October 1996 - June 2002                       *
!*                                                                     *
!*    This program contains the subroutines that advect and diffuse    *
!*  tracers horizontally, along with subroutines that handle the       *
!*  registration of tracers and related subroutines.                   *
!*                                                                     *
!*    advect_tracer advects tracer concentrations using a combination  *
!*  of the modified flux advection scheme from Easter (Mon. Wea. Rev., *
!*  1993) with tracer distributions given by the monotonic modified    *
!*  van Leer scheme proposed by Lin et al. (Mon. Wea. Rev., 1994).     *
!*  This scheme conserves the total amount of tracer while avoiding    *
!*  spurious maxima and minima of the tracer concentration.  If a      *
!*  higher order accuracy scheme is needed, I would suggest the mono-  *
!*  tonic piecewise parabolic method, as described in Carpenter et al. *
!*  (MWR, 1990).  advect_tracer has 4 arguments, described below. This *
!*  subroutine determines the volume of a layer in a grid cell at the  *
!*  previous instance when the tracer concentration was changed, so    *
!*  it is essential that the volume fluxes should be correct.  It is   *
!*  also important that the tracer advection occurs before each        *
!*  calculation of the diabatic forcing.                               *
!*                                                                     *
!*    In addition, each of the tracers are subject to Fickian along-   *
!*  coordinate diffusion if Khtr is defined and positive.  The tracer  *
!*  diffusion uses a suitable number of iterations to guarantee        *
!*  stability with an arbitrarily large time step.  tracer_hordiff is  *
!*  called by advect_tracer.                                           *
!*                                                                     *
!*    This file also contains register_tracer, which is called to      *
!*  indicate the tracers that will be advected and diffused.           *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, vh                                    *
!*    j    x ^ x ^ x   At >:  u, uh                                    *
!*    j    > o > o >   At o:  tr, h                                    *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ptrs
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_domains, only : pass_var, pass_vector, sum_across_PEs, max_across_PEs
use MOM_checksums, only : hchksum
use MOM_EOS, only : calculate_density
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE_types, only : MEKE_type
use MOM_variables, only : ocean_OBC_type, thermo_var_ptrs, OBC_FLATHER_E
use MOM_variables, only : OBC_FLATHER_W, OBC_FLATHER_N, OBC_FLATHER_S

implicit none ; private

#include <MOM_memory.h>

public advect_tracer, tracer_hordiff, register_tracer, advect_tracer_init
public add_tracer_diagnostics, add_tracer_2d_diagnostics, add_tracer_OBC_values
public tracer_vertdiff, advect_tracer_diag_init, tracer_end

type :: tracer
  real, dimension(:,:,:), pointer :: t => NULL()
                     ! The array containing the tracer concentration.
  real :: OBC_inflow_conc = 0.0  ! A tracer concentration for generic inflows.
  real, dimension(:,:,:), pointer :: OBC_in_u => NULL(), OBC_in_v => NULL()
             ! These arrays contain structured values for flow into the domain
             ! that are specified in open boundary conditions through u- and
             ! v- faces of the tracer cell.
  real, dimension(:,:,:), pointer :: ad_x => NULL(), ad_y => NULL()
             ! The arrays in which x- & y- advective fluxes are stored.
  real, dimension(:,:,:), pointer :: df_x => NULL(), df_y => NULL()
             ! The arrays in which x- & y- diffusive fluxes are stored.
  real, dimension(:,:), pointer :: ad2d_x => NULL(), ad2d_y => NULL()
             ! The arrays in which vertically summed x- & y- advective fluxes
             ! are stored in units of CONC m3 s-1..
  real, dimension(:,:), pointer :: df2d_x => NULL(), df2d_y => NULL()
             ! The arrays in which vertically summed x- & y- diffusive fluxes
             ! are stored in units of CONC m3 s-1..
  character(len=32) :: name  ! A tracer name for error messages.
end type tracer

type, public :: advect_tracer_CS ; private
  real    :: dt             ! The baroclinic dynamics time step, in s.
  real    :: KhTr           ! The along-isopycnal tracer diffusivity in m2/s.
  real    :: KhTr_Slope_Cff ! The non-dimensional coefficient in KhTr formula
  real    :: KhTr_min       ! Minimum along-isopycnal tracer diffusivity in m2/s.
  real    :: KhTr_max       ! Maximum along-isopycnal tracer diffusivity in m2/s.
  real    :: KhTr_passivity_coeff ! Passivity coefficient that scales Rd/dx (default = 0)
                            ! where passivity is the ratio between along-isopycnal
                            ! tracer mixing and thickness mixing
  real    :: KhTr_passivity_min   ! Passivity minimum (default = 1/2)
  real    :: ML_KhTR_scale  ! With Diffuse_ML_interior, the ratio of the truly
                            ! horizontal diffusivity in the mixed layer to the
                            ! epipycnal diffusivity.  Nondim.
  logical :: Diffuse_ML_interior  ! If true, diffuse along isopycnals between
                            ! the mixed layer and the interior.
  logical :: check_diffusive_CFL  ! If true, automatically iterate the diffusion
                            ! to ensure that the diffusive equivalent of the CFL
                            ! limit is not violated.
  integer :: ntr = 0        ! The number of registered tracers.
  type(tracer) :: Tr(MAX_FIELDS)  ! The array of registered tracers.
  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                            ! ocean diagnostic fields and control variables.
  logical :: debug           ! If true, write verbose checksums for debugging purposes.
  integer :: id_KhTr_u, id_KhTr_v
end type advect_tracer_CS

type p2d
  real, dimension(:,:), pointer :: p => NULL()
end type p2d
type p2di
  integer, dimension(:,:), pointer :: p => NULL()
end type p2di

integer :: id_clock_advect, id_clock_diffuse, id_clock_epimix
integer :: id_clock_pass, id_clock_sync

contains

subroutine register_tracer(tr1, name, param_file, CS, ad_x, ad_y, &
                           df_x, df_y, OBC_inflow, OBC_in_u, OBC_in_v)
  real, dimension(NXMEM_,NYMEM_,NKMEM_), target :: tr1
  character(len=*), intent(in)               :: name
  type(param_file_type), intent(in)          :: param_file
  type(advect_tracer_CS), pointer            :: CS
  real, pointer, dimension(:,:,:), optional  :: ad_x, ad_y, df_x, df_y
  real, intent(in), optional                 :: OBC_inflow
  real, pointer, dimension(:,:,:), optional  :: OBC_in_u, OBC_in_v
! This subroutine registers a tracer to be advected and horizontally
! diffused.

! Arguments: tr1 - The pointer to the tracer, in arbitrary concentration units (CONC).
!  (in)      name - The name to be used in messages about the tracer.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      ad_x - An array in which zonal advective fluxes are stored in
!                   units of CONC m3 s-1.
!  (in)      ad_y - An array in which meridional advective fluxes are stored
!                   in units of CONC m3 s-1.
!  (in)      df_x - An array in which zonal diffusive fluxes are stored in
!                   units of CONC m3 s-1.
!  (in)      df_y - An array in which meridional diffusive fluxes are stored
!                   in units of CONC m3 s-1.
!  (in)      OBC_inflow - The value of the tracer for all inflows via the open
!                         boundary conditions for which OBC_in_u or OBC_in_v are
!                         not specified, in the same units as tr (CONC).
!  (in)      OBC_in_u - The value of the tracer at inflows through u-faces of
!                       tracer cells, in the same units as tr (CONC).
!  (in)      OBC_in_v - The value of the tracer at inflows through v-faces of
!                       tracer cells, in the same units as tr (CONC).

  integer :: m, ntr
  type(tracer) :: temp
  character(len=256) :: mesg    ! Message for error messages.

  if (.not. associated(CS)) call advect_tracer_init(param_file, CS)

  if (CS%ntr>=MAX_FIELDS) then
    write(mesg,'("Increase MAX_FIELDS in MOM_memory.h to at least ",I3," to allow for &
        &all the tracers being registered via register_tracer.")') CS%ntr+1
    call MOM_error(FATAL,"MOM register_tracer: "//mesg)
  endif
  CS%ntr = CS%ntr + 1
  ntr = CS%ntr

  CS%Tr(ntr)%name = trim(name)
  CS%Tr(ntr)%t => tr1

  if (present(ad_x)) then ; if (associated(ad_x)) CS%Tr(ntr)%ad_x => ad_x ; endif
  if (present(ad_y)) then ; if (associated(ad_y)) CS%Tr(ntr)%ad_y => ad_y ; endif
  if (present(df_x)) then ; if (associated(df_x)) CS%Tr(ntr)%df_x => df_x ; endif
  if (present(df_y)) then ; if (associated(df_y)) CS%Tr(ntr)%df_y => df_y ; endif
  if (present(OBC_inflow)) CS%Tr(ntr)%OBC_inflow_conc = OBC_inflow
  if (present(OBC_in_u)) then ; if (associated(OBC_in_u)) &
                                    CS%Tr(ntr)%OBC_in_u => OBC_in_u ; endif
  if (present(OBC_in_v)) then ; if (associated(OBC_in_v)) &
                                    CS%Tr(ntr)%OBC_in_v => OBC_in_v ; endif

end subroutine register_tracer

subroutine add_tracer_OBC_values(name, CS, OBC_inflow, OBC_in_u, OBC_in_v)
  character(len=*), intent(in)               :: name
  type(advect_tracer_CS), pointer            :: CS
  real, intent(in), optional                 :: OBC_inflow
  real, pointer, dimension(:,:,:), optional  :: OBC_in_u, OBC_in_v
! This subroutine adds open boundary condition concentrations for a tracer that
! has previously been registered by a call to register_tracer.

! Arguments: name - The name of the tracer for which the diagnostic pointers.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      OBC_inflow - The value of the tracer for all inflows via the open
!                         boundary conditions for which OBC_in_u or OBC_in_v are
!                         not specified, in the same units as tr (CONC).
!  (in)      OBC_in_u - The value of the tracer at inflows through u-faces of
!                       tracer cells, in the same units as tr (CONC).
!  (in)      OBC_in_v - The value of the tracer at inflows through v-faces of
  integer :: m

  if (.not. associated(CS)) call MOM_error(FATAL, "add_tracer_OBC_values :"// &
       "register_tracer must be called before add_tracer_OBC_values")

  do m=1,CS%ntr ; if (CS%Tr(m)%name == trim(name)) exit ; enddo

  if (m <= CS%ntr) then
    if (present(OBC_inflow)) CS%Tr(m)%OBC_inflow_conc = OBC_inflow
    if (present(OBC_in_u)) then ; if (associated(OBC_in_u)) &
                                      CS%Tr(m)%OBC_in_u => OBC_in_u ; endif
    if (present(OBC_in_v)) then ; if (associated(OBC_in_v)) &
                                      CS%Tr(m)%OBC_in_v => OBC_in_v ; endif
  else
    call MOM_error(FATAL, "MOM_tracer: register_tracer must be called for "//&
             trim(name)//" before add_tracer_OBC_values is called for it.")
  endif

end subroutine add_tracer_OBC_values

subroutine add_tracer_diagnostics(name, CS, ad_x, ad_y, df_x, df_y)
  character(len=*), intent(in)   :: name
  type(advect_tracer_CS), pointer :: CS
  real, pointer, dimension(:,:,:), optional :: ad_x, ad_y, df_x, df_y
! This subroutine adds diagnostic arrays for a tracer that has previously been
! registered by a call to register_tracer.

! Arguments: name - The name of the tracer for which the diagnostic pointers.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      ad_x - An array in which zonal advective fluxes are stored in
!                   units of CONC m3 s-1.
!  (in)      ad_y - An array in which meridional advective fluxes are stored
!                   in units of CONC m3 s-1.
!  (in)      df_x - An array in which zonal diffusive fluxes are stored in
!                   units of CONC m3 s-1.
!  (in)      df_y - An array in which meridional diffusive fluxes are stored
!                   in units of CONC m3 s-1.
  integer :: m
  logical :: write_warning

  if (.not. associated(CS)) call MOM_error(FATAL, "add_tracer_diagnostics: "// &
       "register_tracer must be called before add_tracer_diagnostics")

  do m=1,CS%ntr ; if (CS%Tr(m)%name == trim(name)) exit ; enddo

  if (m <= CS%ntr) then
    if (present(ad_x)) then ; if (associated(ad_x)) CS%Tr(m)%ad_x => ad_x ; endif
    if (present(ad_y)) then ; if (associated(ad_y)) CS%Tr(m)%ad_y => ad_y ; endif
    if (present(df_x)) then ; if (associated(df_x)) CS%Tr(m)%df_x => df_x ; endif
    if (present(df_y)) then ; if (associated(df_y)) CS%Tr(m)%df_y => df_y ; endif
  else
    call MOM_error(FATAL, "MOM_tracer: register_tracer must be called for "//&
             trim(name)//" before add_tracer_diagnostics is called for it.")
  endif
  if (CS%Diffuse_ML_interior) then
    write_warning = .false.
    if (present(df_x)) then ; if (associated(df_x)) write_warning = .true. ; endif
    if (present(df_y)) then ; if (associated(df_y)) write_warning = .true. ; endif
    if (write_warning .and. is_root_pe()) call MOM_error(WARNING, &
      "add_tracer_diagnostics: "// &
      "3-d tracer diffusion diagnostics are not complete with DIFFUSE_ML_TO_INTERIOR "// &
      "defined.  Use 2-d tracer diffusion diagnostics instead for total fluxes.")
  endif 

end subroutine add_tracer_diagnostics

subroutine add_tracer_2d_diagnostics(name, CS, ad_x, ad_y, df_x, df_y)
  character(len=*), intent(in)   :: name
  type(advect_tracer_CS), pointer :: CS
  real, pointer, dimension(:,:), optional :: ad_x, ad_y, df_x, df_y
! This subroutine adds diagnostic arrays for a tracer that has previously been
! registered by a call to register_tracer.

! Arguments: name - The name of the tracer for which the diagnostic pointers.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      ad_x - An array in which the vertically summed zonal advective
!                   fluxes are stored in units of CONC m3 s-1.
!  (in)      ad_y - An array in which the vertically summed meridional advective
!                   fluxes are stored in units of CONC m3 s-1.
!  (in)      df_x - An array in which the vertically summed zonal diffusive
!                   fluxes are stored in units of CONC m3 s-1.
!  (in)      df_y - An array in which the vertically summed meridional diffusive
!                   fluxes are stored in units of CONC m3 s-1.
  integer :: m

  if (.not. associated(CS)) call MOM_error(FATAL, "add_tracer_diagnostics: "// &
       "register_tracer must be called before add_tracer_diagnostics")

  do m=1,CS%ntr ; if (CS%Tr(m)%name == trim(name)) exit ; enddo

  if (m <= CS%ntr) then
    if (present(ad_x)) then ; if (associated(ad_x)) CS%Tr(m)%ad2d_x => ad_x ; endif
    if (present(ad_y)) then ; if (associated(ad_y)) CS%Tr(m)%ad2d_y => ad_y ; endif
    if (present(df_x)) then ; if (associated(df_x)) CS%Tr(m)%df2d_x => df_x ; endif
    if (present(df_y)) then ; if (associated(df_y)) CS%Tr(m)%df2d_y => df_y ; endif
  else
    call MOM_error(FATAL, "MOM_tracer: register_tracer must be called for "//&
             trim(name)//" before add_tracer_2d_diagnostics is called for it.")
  endif

end subroutine add_tracer_2d_diagnostics


subroutine advect_tracer(h_end, uhtr, vhtr, OBC, dt, G, CS)
  real, dimension(NXMEM_,NYMEM_,NKMEM_),  intent(in)    :: h_end
  real, dimension(NXMEMQ_,NYMEM_,NKMEM_), intent(in)    :: uhtr
  real, dimension(NXMEM_,NYMEMQ_,NKMEM_), intent(in)    :: vhtr
  type(ocean_OBC_type),                   pointer       :: OBC
  real,                                   intent(in)    :: dt
  type(ocean_grid_type),                  intent(inout) :: G
  type(advect_tracer_CS),                 pointer       :: CS
!    This subroutine time steps the tracer concentration.
!  A monotonic, conservative, weakly diffusive scheme is used.

! Arguments: h_end - Layer thickness after advection, in m or kg m-2.
!  (in)      uhtr - Accumulated volume or mass fluxes through zonal faces,
!                   in m3 or kg.
!  (in)      vhtr - Accumulated volume or mass fluxes through meridional faces,
!                   in m3 or kg.
!  (in)      OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 advect_tracer_init.

  type(tracer) :: Tr(MAX_FIELDS) ! The array of registered tracers.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    hprev           ! The cell volume at the end of the previous tracer
                    ! change, in m3.
  real, dimension(SZIQ_(G),SZJ_(G),SZK_(G)) :: &
    uhr             ! The remaining zonal thickness flux, in m3.
  real, dimension(SZI_(G),SZJQ_(G),SZK_(G)) :: &
    vhr             ! The remaining meridional thickness fluxes, in m3.
  real :: uh_neglect(SZIQ_(G),SZJ_(G)) ! uh_neglect and vh_neglect are the
  real :: vh_neglect(SZI_(G),SZJQ_(G)) ! magnitude of remaining transports that
                                ! can be simply discarded, in m3 or kg.

  real :: landvolfill         ! An arbitrary? nonzero cell volume, m3.
  real :: Idt                 ! 1/dt in s-1.
  logical :: domore_u(SZJ_(G),SZK_(G))  ! domore__ indicate whether there is more
  logical :: domore_v(SZJQ_(G),SZK_(G)) ! advection to be done in the corresponding
                                ! row or column.
  logical :: x_first            ! If true, advect in the x-direction first.
  integer :: max_iter           ! The maximum number of iterations in
                                ! each layer.
  integer :: domore_k(SZK_(G))
  integer :: stensil            ! The stensil of the advection scheme.
  integer :: nsten_halo         ! The number of stensils that fit in the halos.
  integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz, itt, ntr, do_any
  integer :: isv, iev, jsv, jev ! The valid range of the indices.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  landvolfill = 1.0e-20         ! This is arbitrary, but must be positive.
  stensil = 2                   ! The scheme's stensil; 2 for PLM.

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_tracer: "// &
       "register_tracer must be called before advect_tracer.")
  if (CS%ntr==0) return
  call cpu_clock_begin(id_clock_advect)
  x_first = (MOD(G%first_direction,2) == 0)

  ntr = CS%ntr
  do m=1,ntr ; Tr(m) = CS%Tr(m) ; enddo
  Idt = 1.0/dt

  max_iter = 2*INT(CEILING(dt/CS%dt)) + 1

! This initializes the halos of uhr and vhr because pass_vector might do
! calculations on them, even though they are never used.
  uhr(:,:,:) = 0.0 ; vhr(:,:,:) = 0.0
  hprev(:,:,:) = landvolfill

  do k=1,nz
    domore_k(k)=1
!  Put the remaining (total) thickness fluxes into uhr and vhr.
    do j=js,je ; do I=is-1,ie ; uhr(I,j,k) = uhtr(I,j,k) ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; vhr(i,J,k) = vhtr(i,J,k) ; enddo ; enddo
!   This loop reconstructs the thickness field the last time that the
! tracers were updated, probably just after the diabatic forcing.  A useful
! diagnostic could be to compare this reconstruction with that older value.
    do i=is,ie ; do j=js,je
      hprev(i,j,k) = max(0.0, G%DXDYh(i,j)*h_end(i,j,k) + &
           ((uhr(I,j,k) - uhr(I-1,j,k)) + (vhr(i,J,k) - vhr(i,J-1,k))))
! In the case that the layer is now dramatically thinner than it was previously,
! add a bit of mass to avoid truncation errors.  This will lead to
! non-conservation of tracers 
      hprev(i,j,k) = hprev(i,j,k) + &
                     max(0.0, 1.0e-13*hprev(i,j,k) - G%DXDYh(i,j)*h_end(i,j,k))
    enddo ; enddo
  enddo
  do j=jsd,jed ; do I=isd,ied-1
    uh_neglect(I,j) = G%H_subroundoff*MIN(G%DXDYh(i,j),G%DXDYh(i+1,j))
  enddo ; enddo
  do J=jsd,jed-1 ; do i=isd,ied
    vh_neglect(i,J) = G%H_subroundoff*MIN(G%DXDYh(i,j),G%DXDYh(i,j+1))
  enddo ; enddo

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
    if (associated(Tr(m)%ad2d_x)) then
      do j=jsd,jed ; do i=isd,ied ; Tr(m)%ad2d_x(I,j) = 0.0 ; enddo ; enddo
    endif
    if (associated(Tr(m)%ad2d_y)) then
      do J=jsd,jed ; do i=isd,ied ; Tr(m)%ad2d_y(i,J) = 0.0 ; enddo ; enddo
    endif
  enddo

  isv = is ; iev = ie ; jsv = js ; jev = je

  do itt=1,max_iter

    if (isv > is-stensil) then
      call cpu_clock_begin(id_clock_pass)
      call pass_vector(uhr, vhr, G%Domain)
      do m=1,ntr ; call pass_var(Tr(m)%t, G%Domain, complete=.false.) ; enddo
      call pass_var(hprev, G%Domain)
      call cpu_clock_end(id_clock_pass)

      nsten_halo = min(is-isd,ied-ie,js-jsd,jed-je)/stensil
      isv = is-nsten_halo*stensil ; jsv = js-nsten_halo*stensil
      iev = ie+nsten_halo*stensil ; jev = je+nsten_halo*stensil
      ! Reevaluate domore_u & domore_v unless the valid range is the same size as
      ! before.  Also, do this if there is Strang splitting.
      if ((nsten_halo > 1) .or. (itt==1)) then
        do k=1,nz ; if (domore_k(k) > 0) then
          do j=jsv,jev ; if (.not.domore_u(j,k)) then
            do i=isv+stensil-1,iev-stensil; if (uhr(I,j,k) /= 0.0) then
              domore_u(j,k) = .true. ; exit
            endif ; enddo ! i-loop
          endif ; enddo
          do J=jsv+stensil-1,jev-stensil ; if (.not.domore_v(J,k)) then
            do i=isv+stensil,iev-stensil; if (vhr(i,J,k) /= 0.0) then
              domore_v(J,k) = .true. ; exit
            endif ; enddo ! i-loop
          endif ; enddo
          
          !   At this point, domore_k is global.  Change it so that it indicates
          ! whether any work is needed on a layer on this processor.
          domore_k(k) = 0
          do j=jsv,jev ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
          do J=jsv+stensil-1,jev-stensil ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo

        endif ; enddo ! k-loop    
      endif
    endif
    
    ! Set the range of valid points after this iteration.
    isv = isv + stensil ; iev = iev - stensil
    jsv = jsv + stensil ; jev = jev - stensil

    do k=1,nz ; if (domore_k(k) > 0) then
!    To ensure positive definiteness of the thickness at each iteration, the
!  mass fluxes out of each layer are checked each step, and limited to keep
!  the thicknesses positive.  This means that several iteration may be required
!  for all the transport to happen.  The sum over domore_k keeps the processors
!  synchronized.  This may not be very efficient, but it should be reliable.

      if (x_first) then
  !    First, advect zonally.
        call advect_x(Tr, hprev, uhr, uh_neglect, OBC, domore_u, ntr, Idt, &
                      isv, iev, jsv-stensil, jev+stensil, k, G)

  !    Next, advect meridionally.
        call advect_y(Tr, hprev, vhr, vh_neglect, OBC, domore_v, ntr, Idt, &
                      isv, iev, jsv, jev, k, G)

        domore_k(k) = 0
        do j=jsv-stensil,jev+stensil ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
        do J=jsv-1,jev ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo
      else
  !    First, advect meridionally.
        call advect_y(Tr, hprev, vhr, vh_neglect, OBC, domore_v, ntr, Idt, &
                      isv-stensil, iev+stensil, jsv, jev, k, G)

  !    Next, advect zonally.
        call advect_x(Tr, hprev, uhr, uh_neglect, OBC, domore_u, ntr, Idt, &
                      isv, iev, jsv, jev, k, G)

        domore_k(k) = 0
        do j=jsv,jev ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
        do J=jsv-1,jev ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo
      endif

    endif ; enddo ! End of k-loop

    ! If the advection just isn't finishing after max_iter, move on.
    if (itt >= max_iter) exit

    ! Exit if there are no layers that need more iterations.
    if (isv > is-stensil) then
      do_any = 0
      call cpu_clock_begin(id_clock_sync)
      call sum_across_PEs(domore_k(:), nz)
      call cpu_clock_end(id_clock_sync)
      do k=1,nz ; do_any = do_any + domore_k(k) ; enddo
      if (do_any == 0) exit
    endif

  enddo ! Iterations loop

  call cpu_clock_end(id_clock_advect)

end subroutine advect_tracer

subroutine advect_x(Tr, hprev, uhr, uh_neglect, OBC, domore_u, ntr, Idt, &
                    is, ie, js, je, k, G)
  type(tracer), dimension(ntr),           intent(inout) :: Tr
  real, dimension(NXMEM_,NYMEM_,NKMEM_),  intent(inout) :: hprev
  real, dimension(NXMEMQ_,NYMEM_,NKMEM_), intent(inout) :: uhr
  real, dimension(NXMEMQ_,NYMEM_),        intent(inout) :: uh_neglect
  type(ocean_OBC_type),                   pointer       :: OBC
  logical, dimension(NYMEM_,NKMEM_),      intent(inout) :: domore_u
  real,                                   intent(in)    :: Idt
  integer,                                intent(in)    :: ntr, is, ie, js, je,k
  type(ocean_grid_type),                  intent(inout) :: G
  !   This subroutine does 1-d flux-form advection in the zonal direction using
  ! a monotonic piecewise linear scheme.
  real, dimension(SZIQ_(G),ntr) :: &
    slope_x, &      ! The concentration slope per grid point in units of
                    ! concentration (nondim.).
    flux_x          ! The tracer flux across a boundary in m3*conc or kg*conc.
  real :: maxslope            ! The maximum concentration slope per grid point
                              ! consistent with monotonicity, in conc. (nondim.).
  real :: hup, hlos           ! hup is the upwind volume, hlos is the
                              ! part of that volume that might be lost
                              ! due to advection out the other side of
                              ! the grid box, both in m3 or kg.
  real :: uhh(SZIQ_(G))       ! The zonal flux that occurs during the
                              ! current iteration, in m3 or kg.
  real, dimension(SZIQ_(G)) :: &
    hlst, Ihnew, &      ! Work variables with units of m3 or kg and m-3 or kg-1.
    ts2                 ! A nondimensional work variable.
  real :: min_h         ! The minimum thickness that can be realized during
                        ! any of the passes, in m or kg m-2.
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected, in m.
  logical :: do_i(SZIQ_(G))     ! If true, work on given points.
  logical :: do_any_i
  integer :: i, j, m

  min_h = 0.1*G%Angstrom
  h_neglect = G%H_subroundoff

  do I=is-1,ie ; ts2(I) = 0.0 ; enddo
  
  do j=js,je ; if (domore_u(j,k)) then
    domore_u(j,k) = .false.

!   Calculate the i-direction profiles (slopes) of each tracer that
! is being advected.
    do m=1,ntr ; do i=is-1,ie+1
      if (ABS(Tr(m)%t(i+1,j,k)-Tr(m)%t(i,j,k)) < &
          ABS(Tr(m)%t(i,j,k)-Tr(m)%t(i-1,j,k))) then
        maxslope = 4.0*(Tr(m)%t(i+1,j,k)-Tr(m)%t(i,j,k))
      else
        maxslope = 4.0*(Tr(m)%t(i,j,k)-Tr(m)%t(i-1,j,k))
      endif
      if ((Tr(m)%t(i+1,j,k)-Tr(m)%t(i,j,k)) * (Tr(m)%t(i,j,k)-Tr(m)%t(i-1,j,k)) < 0.0) then
        slope_x(i,m) = 0.0
      elseif (ABS(Tr(m)%t(i+1,j,k)-Tr(m)%t(i-1,j,k))<ABS(maxslope)) then
        slope_x(i,m) = G%umask(I,j)*G%umask(I-1,j) * &
                       0.5*(Tr(m)%t(i+1,j,k)-Tr(m)%t(i-1,j,k))
      else
        slope_x(i,m) = G%umask(I,j)*G%umask(I-1,j) * 0.5*maxslope
      endif
    enddo ; enddo

!   Calculate the i-direction fluxes of each tracer, using as much
! the minimum of the remaining mass flux (uhr) and the half the mass
! in the cell plus whatever part of its half of the mass flux that
! the flux through the other side does not require.
    do I=is-1,ie
      if (uhr(I,j,k) == 0.0) then
        uhh(I) = 0.0
      elseif (uhr(I,j,k) < 0.0) then
        hup = (hprev(i+1,j,k)-G%DXDYh(i+1,j)*G%Angstrom*0.1)
        hlos = MAX(0.0,uhr(I+1,j,k))
        if (((hup + uhr(I,j,k) - hlos) < 0.0) .and. &
            ((0.5*hup + uhr(I,j,k)) < 0.0)) then
          uhh(I) = MIN(-0.5*hup,-hup+hlos,0.0)
          domore_u(j,k) = .true.
        else
          uhh(I) = uhr(I,j,k)
        endif
        ts2(I) = 0.5*(1.0 + uhh(I)/(hprev(i+1,j,k)+h_neglect))
      else
        hup = (hprev(i,j,k)-G%DXDYh(i,j)*G%Angstrom*0.1)
        hlos = MAX(0.0,-uhr(I-1,j,k))
        if (((hup - uhr(I,j,k) - hlos) < 0.0) .and. &
            ((0.5*hup - uhr(I,j,k)) < 0.0)) then
          uhh(I) = MAX(0.5*hup,hup-hlos,0.0)
          domore_u(j,k) = .true.
        else
          uhh(I) = uhr(I,j,k)
        endif
        ts2(I) = 0.5*(1.0 - uhh(I)/(hprev(i,j,k)+h_neglect))
      endif
    enddo
    do m=1,ntr ; do I=is-1,ie
      if (uhh(I) >= 0.0) then
        flux_x(I,m) = uhh(I)*(Tr(m)%t(i,j,k) + slope_x(i,m)*ts2(I))
      else
        flux_x(I,m) = uhh(I)*(Tr(m)%t(i+1,j,k) - slope_x(i+1,m)*ts2(I))
      endif
    enddo ; enddo
    if (associated(OBC)) then ; if (OBC%apply_OBC_u) then
      do_any_i = .false.
      do I=is-1,ie
        do_i(I) = .false.
        if (OBC%OBC_mask_u(I,j) .and. uhr(I,j,k) /= 0.0) then
          ! Tracer fluxes are set to prescribed values only for inflows
          ! from masked areas.
          if (((uhr(I,j,k) > 0.0) .and. ((G%hmask(i,j) < 0.5) .or. &
                  (OBC%OBC_kind_u(I,j) == OBC_FLATHER_W))) .or. &
              ((uhr(I,j,k) < 0.0) .and. ((G%hmask(i+1,j) < 0.5) .or. &
                  (OBC%OBC_kind_u(I,j) == OBC_FLATHER_E))) ) then
            do_i(I) = .true. ; do_any_i = .true.
            uhh(I) = uhr(I,j,k)
          endif
        endif
      enddo
      if (do_any_i) then ; do m=1,ntr ; do i=is,ie ; if (do_i(i)) then
        if (associated(Tr(m)%OBC_in_u)) then
          flux_x(I,m) = uhh(I)*Tr(m)%OBC_in_u(I,j,k)
        else ; flux_x(I,m) = uhh(I)*Tr(m)%OBC_inflow_conc ; endif
      endif ; enddo ; enddo ; endif
    endif ; endif

!   Calculate new tracer concentration in each cell after accounting
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
        elseif (hprev(i,j,k) < h_neglect*G%DXDYh(i,j)) then
          hlst(i) = hlst(i) + (h_neglect*G%DXDYh(i,j) - hprev(i,j,k))
          Ihnew(i) = 1.0 / (h_neglect*G%DXDYh(i,j))
        else ;  Ihnew(i) = 1.0 / hprev(i,j,k) ; endif
      else
        do_i(i) = .false.
      endif
    enddo
    do m=1,ntr
      do i=is,ie ; if ((do_i(i)) .and. (Ihnew(i) > 0.0)) then
        Tr(m)%t(i,j,k) = (Tr(m)%t(i,j,k) * hlst(i) - &
                          (flux_x(I,m) - flux_x(I-1,m))) * Ihnew(i)
      endif ; enddo
      if (associated(Tr(m)%ad_x)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad_x(I,j,k) = Tr(m)%ad_x(I,j,k) + flux_x(I,m)*Idt
      endif ; enddo ; endif
      if (associated(Tr(m)%ad2d_x)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad2d_x(I,j) = Tr(m)%ad2d_x(I,j) + flux_x(I,m)*Idt
      endif ; enddo ; endif
    enddo

  endif ; enddo ! End of j-loop.

end subroutine advect_x

subroutine advect_y(Tr, hprev, vhr, vh_neglect, OBC, domore_v, ntr, Idt, &
                    is, ie, js, je, k, G)
  type(tracer), dimension(ntr),           intent(inout) :: Tr
  real, dimension(NXMEM_,NYMEM_,NKMEM_),  intent(inout) :: hprev
  real, dimension(NXMEM_,NYMEMQ_,NKMEM_), intent(inout) :: vhr
  real, dimension(NXMEM_,NYMEMQ_),        intent(inout) :: vh_neglect
  type(ocean_OBC_type),                   pointer       :: OBC
  logical, dimension(NYMEMQ_,NKMEM_),     intent(inout) :: domore_v
  real,                                   intent(in)    :: Idt
  integer,                                intent(in)    :: ntr, is, ie, js, je,k
  type(ocean_grid_type),                  intent(inout) :: G
  !   This subroutine does 1-d flux-form advection using a monotonic piecewise
  ! linear scheme.
  real, dimension(SZI_(G),ntr,SZJQ_(G)) :: &
    slope_y, &      ! The concentration slope per grid point in units of
                    ! concentration (nondim.).
    flux_y          ! The tracer flux across a boundary in m3 * conc or kg*conc.
  real :: maxslope            ! The maximum concentration slope per grid point
                              ! consistent with monotonicity, in conc. (nondim.).
  real :: vhh(SZI_(G),SZJQ_(G)) ! The meridional flux that occurs during the
                              ! current iteration, in m3 or kg.
  real :: hup, hlos           ! hup is the upwind volume, hlos is the
                              ! part of that volume that might be lost
                              ! due to advection out the other side of
                              ! the grid box, both in m3 or kg.
  real, dimension(SZIQ_(G)) :: &
    hlst, Ihnew, &      ! Work variables with units of m3 or kg and m-3 or kg-1.
    ts2                 ! A nondimensional work variable.
  real :: min_h         ! The minimum thickness that can be realized during
                        ! any of the passes, in m or kg m-2.
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected, in m.
  logical :: do_j_tr(SZJ_(G))   ! If true, calculate the tracer profiles.
  logical :: do_i(SZIQ_(G))     ! If true, work on given points.
  logical :: do_any_i
  integer :: i, j, m

  min_h = 0.1*G%Angstrom
  h_neglect = G%H_subroundoff

  do i=is,ie ; ts2(i) = 0.0 ; enddo
!   Calculate the j-direction profiles (slopes) of each tracer that
! is being advected.
  do_j_tr(js-1) = domore_v(js-1,k) ; do_j_tr(je+1) = domore_v(je,k)
  do j=js,je ; do_j_tr(j) = (domore_v(J-1,k) .or. domore_v(J,k)) ; enddo
  do j=js-1,je+1 ; if (do_j_tr(j)) then ; do m=1,ntr ; do i=is,ie
    if (ABS(Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j,k)) < &
        ABS(Tr(m)%t(i,j,k)-Tr(m)%t(i,j-1,k))) then
      maxslope = 4.0*(Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j,k))
    else
      maxslope = 4.0*(Tr(m)%t(i,j,k)-Tr(m)%t(i,j-1,k))
    endif

    if ((Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j,k))*(Tr(m)%t(i,j,k)-Tr(m)%t(i,j-1,k)) < 0.0) then
      slope_y(i,m,j) = 0.0
    elseif (ABS(Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j-1,k))<ABS(maxslope)) then
      slope_y(i,m,j) = G%vmask(i,J) * G%vmask(i,J-1) * &
                     0.5*(Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j-1,k))
    else
      slope_y(i,m,j) = G%vmask(i,J) * G%vmask(i,J-1) * 0.5*maxslope
    endif
  enddo ; enddo ; endif ; enddo ! End of i-, m-, & j- loops.

!   Calculate the j-direction fluxes of each tracer, using as much
! the minimum of the remaining mass flux (vhr) and the half the mass
! in the cell plus whatever part of its half of the mass flux that
! the flux through the other side does not require.
  do J=js-1,je ; if (domore_v(J,k)) then
    domore_v(J,k) = .false.
    do i=is,ie
      if (vhr(i,J,k) == 0.0) then
        vhh(i,J) = 0.0
      elseif (vhr(i,J,k) < 0.0) then
        hup = (hprev(i,j+1,k)-G%DXDYh(i,j+1)*G%Angstrom*0.1)
        hlos = MAX(0.0,vhr(i,J+1,k))
        if ((((hup - hlos) + vhr(i,J,k)) < 0.0) .and. &
            ((0.5*hup + vhr(i,J,k)) < 0.0)) then
          vhh(i,J) = MIN(-0.5*hup,-hup+hlos,0.0)
          domore_v(J,k) = .true.
        else
          vhh(i,J) = vhr(i,J,k)
        endif
        ts2(i) = 0.5*(1.0 + vhh(i,J) / (hprev(i,j+1,k)+h_neglect))
      else
        hup = (hprev(i,j,k)-G%DXDYh(i,j)*G%Angstrom*0.1)
        hlos = MAX(0.0,-vhr(i,J-1,k))
        if ((((hup - hlos) - vhr(i,J,k)) < 0.0) .and. &
            ((0.5*hup - vhr(i,J,k)) < 0.0)) then
          vhh(i,J) = MAX(0.5*hup,hup-hlos,0.0)
          domore_v(J,k) = .true.
        else
          vhh(i,J) = vhr(i,J,k)
        endif
        ts2(i) = 0.5*(1.0 - vhh(i,J) / (hprev(i,j,k)+h_neglect))
      endif
    enddo
    do m=1,ntr ; do i=is,ie
      if (vhh(i,J) >= 0.0) then
        flux_y(i,m,J) = vhh(i,J)*(Tr(m)%t(i,j,k) + slope_y(i,m,j)*ts2(i))
      else
        flux_y(i,m,J) = vhh(i,J)*(Tr(m)%t(i,j+1,k) - slope_y(i,m,j+1)*ts2(i))
      endif
    enddo ; enddo

    if (associated(OBC)) then ; if (OBC%apply_OBC_v) then
      do_any_i = .false.
      do i=is,ie
        do_i(i) = .false.
        if (OBC%OBC_mask_v(i,J) .and. vhr(i,J,k) /= 0.0) then
        ! Tracer fluxes are set to prescribed values only for inflows
        ! from masked areas.
          if (((vhr(i,J,k) > 0.0) .and. ((G%hmask(i,j) < 0.5) .or. &
                  (OBC%OBC_kind_v(i,J) == OBC_FLATHER_S))) .or. &
              ((vhr(i,J,k) < 0.0) .and. ((G%hmask(i,j+1) < 0.5) .or. &
                  (OBC%OBC_kind_v(i,J) == OBC_FLATHER_N))) ) then
            do_i(i) = .true. ; do_any_i = .true.
            vhh(i,J) = vhr(i,J,k)
          endif
        endif
      enddo
      if (do_any_i) then ; do m=1,ntr ; do i=is,ie ; if (do_i(i)) then
        if (associated(Tr(m)%OBC_in_v)) then
          flux_y(i,m,J) = vhh(i,J)*Tr(m)%OBC_in_v(i,J,k)
        else ; flux_y(i,m,J) = vhh(i,J)*Tr(m)%OBC_inflow_conc ; endif
      endif ; enddo ; enddo ; endif
    endif ; endif
  else ! not domore_v.
    do i=is,ie ; vhh(i,J) = 0.0 ; enddo
    do m=1,ntr ; do i=is,ie ; flux_y(i,m,J) = 0.0 ; enddo ; enddo
  endif ; enddo ! End of j-loop

  do J=js-1,je ; do i=is,ie
    vhr(i,J,k) = vhr(i,J,k) - vhh(i,J)
    if (abs(vhr(i,J,k)) < vh_neglect(i,J)) vhr(i,J,k) = 0.0
  enddo ; enddo

!   Calculate new tracer concentration in each cell after accounting
! for the j-direction fluxes.
  do j=js,je ; if (do_j_tr(j)) then
    do i=is,ie
      if ((vhh(i,J) /= 0.0) .or. (vhh(i,J-1) /= 0.0)) then
        do_i(i) = .true.
        hlst(i) = hprev(i,j,k)
        hprev(i,j,k) = max(hprev(i,j,k) - (vhh(i,J) - vhh(i,J-1)), 0.0)
        if (hprev(i,j,k) <= 0.0) then ; do_i(i) = .false.
        elseif (hprev(i,j,k) < h_neglect*G%DXDYh(i,j)) then
          hlst(i) = hlst(i) + (h_neglect*G%DXDYh(i,j) - hprev(i,j,k))
          Ihnew(i) = 1.0 / (h_neglect*G%DXDYh(i,j))
        else ;  Ihnew(i) = 1.0 / hprev(i,j,k) ; endif
      else ; do_i(i) = .false. ; endif
    enddo
    do m=1,ntr
      do i=is,ie ; if (do_i(i)) then
        Tr(m)%t(i,j,k) = (Tr(m)%t(i,j,k) * hlst(i) - &
                          (flux_y(i,m,J) - flux_y(i,m,J-1))) * Ihnew(i)
      endif ; enddo
      if (associated(Tr(m)%ad_y)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad_y(i,J,k) = Tr(m)%ad_y(i,J,k) + flux_y(i,m,J)*Idt
      endif ; enddo ; endif
      if (associated(Tr(m)%ad2d_y)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad2d_y(i,J) = Tr(m)%ad2d_y(i,J) + flux_y(i,m,J)*Idt
      endif ; enddo ; endif
    enddo
  endif ; enddo ! End of j-loop.

end subroutine advect_y

subroutine tracer_hordiff(h, dt, MEKE, VarMix, G, CS, tv)
  real, dimension(NXMEM_,NYMEM_,NKMEM_), intent(in)    :: h
  real,                               intent(in)    :: dt
  type(MEKE_type),                    pointer       :: MEKE
  type(VarMix_CS),                    pointer       :: VarMix
  type(ocean_grid_type),              intent(inout) :: G
  type(advect_tracer_CS),             pointer       :: CS
  type(thermo_var_ptrs),              intent(in)    :: tv

!   This subroutine does along-coordinate diffusion of all tracers,
! using the diffusivity in CS%KhTr.  Multiple iterations are
! used (if necessary) so that there is no limit on the acceptable
! time increment.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (in)      dt - Time increment in s.
!  (in/out)  Tr - An array of all of the registered tracers.
!  (in)      VarMix - A structure with information about horizontal diffusivities.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 advect_tracer_init.
!  (in,opt)  tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs,
!                 and these may (probably will) point to some of the same arrays
!                 as Tr does.  tv is required for epipycnal mixing between the
!                 mixed layer and the interior.
  type(tracer) :: Tr(MAX_FIELDS) ! The array of registered tracers.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    Ihdxdy, &     ! The inverse of the volume or mass of fluid in a layer in a
                  ! grid cell, in m-3 or kg-1.
    Kh_h, &       ! The tracer diffusivity averaged to tracer points, in m2 s-1.
    CFL, &        ! A diffusive CFL number for each cell, nondim.
    dTr           ! The change in a tracer's concentration, in units of
                  ! concentration.
  real, dimension(SZIQ_(G),SZJ_(G)) :: &
    khdt_x, &     ! The value of Khtr*dt times the open face width divided by
                  ! the distance between adjacent tracer points, in m2.
    Coef_x, &     ! The coefficients relating zonal tracer differences
                  ! to time-integrated fluxes, in m3 or kg.
    Kh_u          ! Tracer mixing coefficient at u-points, in m2 s-1.
  real, dimension(SZI_(G),SZJQ_(G)) :: &
    khdt_y, &     ! The value of Khtr*dt times the open face width divided by
                  ! the distance between adjacent tracer points, in m2.
    Coef_y, &     ! The coefficients relating meridional tracer differences
                  ! to time-integrated fluxes, in m3 or kg.
    Kh_v          ! Tracer mixing coefficient at u-points, in m2 s-1.

  real :: max_CFL ! The global maximum of the diffusive CFL number.
  logical :: use_VarMix, Resoln_scaled
  integer :: i, j, k, m, is, ie, js, je, nz, ntr, itt, num_itts
  real :: I_numitts  ! The inverse of the number of iterations, num_itts.
  real :: scale      ! The fraction of khdt_x or khdt_y that is applied in this
                     ! layer for this iteration, nondim.
  real :: Idt        ! The inverse of the time step, in s-1.
  real :: h_neglect  ! A thickness that is so small it is usually lost
                     ! in roundoff and can be neglected, in m.
  real :: Kh_loc     ! The local value of Kh, in m2 s-1.
  real :: Res_Fn     ! The local value of the resolution function, nondim.
  real :: Rd_dx      ! The local value of deformation radius over grid-spacing, nondim.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_tracer: "// &
       "register_tracer must be called before tracer_hordiff.")
  if ((CS%ntr==0) .or. ((CS%KhTr <= 0.0) .and. .not.associated(VarMix)) ) return

  call cpu_clock_begin(id_clock_diffuse)

  ntr = CS%ntr
  do m=1,ntr ; Tr(m) = CS%Tr(m) ; enddo
  Idt = 1.0/dt
  h_neglect = G%H_subroundoff

  if (CS%debug) call MOM_tracer_chksum("Before tracer diffusion ", Tr, ntr, G)

  use_VarMix = .false. ; Resoln_scaled = .false.
  if (Associated(VarMix)) then
    use_VarMix = VarMix%use_variable_mixing
    Resoln_scaled = VarMix%Resoln_scaled_KhTr
  endif

  if (use_VarMix) then
    do j=js,je ; do I=is-1,ie
      Kh_loc = CS%KhTr + CS%KhTr_Slope_Cff*VarMix%L2u(I,j)*VarMix%SN_u(I,j)
      if (associated(MEKE%Kh)) &
        Kh_Loc = Kh_Loc + MEKE%KhTr_fac*sqrt(MEKE%Kh(i,j)*MEKE%Kh(i+1,j))
      if (CS%KhTr_max > 0.) Kh_loc = min(Kh_loc, CS%KhTr_max)
      if (Resoln_scaled) &
        Kh_Loc = Kh_Loc * 0.5*(VarMix%Res_fn_h(i,j) + VarMix%Res_fn_h(i+1,j))
      Kh_u(I,j) = max(Kh_loc, CS%KhTr_min)
      if (CS%KhTr_passivity_coeff>0.) then ! Apply passivity
        Rd_dx=0.5*( VarMix%Rd_dx_h(i,j)+VarMix%Rd_dx_h(i+1,j) ) ! Rd/dx at u-points
        Kh_loc=Kh_u(I,j)*max( CS%KhTr_passivity_min, CS%KhTr_passivity_coeff*Rd_dx )
        if (CS%KhTr_max > 0.) Kh_loc = min(Kh_loc, CS%KhTr_max) ! Re-apply max
        Kh_u(I,j) = max(Kh_loc, CS%KhTr_min) ! Re-apply min
      endif
    enddo ; enddo
    do J=js-1,je ;  do i=is,ie
      Kh_loc = CS%KhTr + CS%KhTr_Slope_Cff*VarMix%L2v(i,J)*VarMix%SN_v(i,J)
      if (associated(MEKE%Kh)) &
        Kh_Loc = Kh_Loc + MEKE%KhTr_fac*sqrt(MEKE%Kh(i,j)*MEKE%Kh(i,j+1))
      if (CS%KhTr_max > 0.) Kh_loc = min(Kh_loc, CS%KhTr_max)
      if (Resoln_scaled) &
        Kh_Loc = Kh_Loc * 0.5*(VarMix%Res_fn_h(i,j) + VarMix%Res_fn_h(i,j+1))
      Kh_v(i,J) = max(Kh_loc, CS%KhTr_min)
      if (CS%KhTr_passivity_coeff>0.) then ! Apply passivity
        Rd_dx=0.5*( VarMix%Rd_dx_h(i,j)+VarMix%Rd_dx_h(i,j+1) ) ! Rd/dx at v-points
        Kh_loc=Kh_v(I,j)*max( CS%KhTr_passivity_min, CS%KhTr_passivity_coeff*Rd_dx )
        if (CS%KhTr_max > 0.) Kh_loc = min(Kh_loc, CS%KhTr_max) ! Re-apply max
        Kh_v(i,J) = max(Kh_loc, CS%KhTr_min) ! Re-apply min
      endif
    enddo ; enddo
    do j=js,je ; do I=is-1,ie
      khdt_x(I,j) = dt*(Kh_u(I,j)*(G%dy_u(I,j)*G%IDXu(I,j)))
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      khdt_y(i,J) = dt*(Kh_v(i,J)*(G%dx_v(i,J)*G%IDYv(i,J)))
    enddo ; enddo
  elseif (Resoln_scaled) then
    do j=js,je ; do I=is-1,ie
      Res_fn = 0.5 * (VarMix%Res_fn_h(i,j) + VarMix%Res_fn_h(i+1,j))
      Kh_u(I,j) = max(CS%KhTr * Res_fn, CS%KhTr_min)
      khdt_x(I,j) = dt*(CS%KhTr*(G%dy_u(I,j)*G%IDXu(I,j))) * Res_fn
    enddo ; enddo
    do J=js-1,je ;  do i=is,ie
      Res_fn = 0.5*(VarMix%Res_fn_h(i,j) + VarMix%Res_fn_h(i,j+1))
      Kh_v(i,J) = max(CS%KhTr * Res_fn, CS%KhTr_min)
      khdt_y(i,J) = dt*(CS%KhTr*(G%dx_v(i,J)*G%IDYv(i,J))) * Res_fn
    enddo ; enddo
  else
    if (CS%id_KhTr_u > 0) then ; do j=js,je ; do I=is-1,ie
      Kh_u(I,j) = CS%KhTr
      khdt_x(I,j) = dt*(CS%KhTr*(G%dy_u(I,j)*G%IDXu(I,j)))
    enddo ; enddo ; else ; do j=js,je ; do I=is-1,ie
      khdt_x(I,j) = dt*(CS%KhTr*(G%dy_u(I,j)*G%IDXu(I,j)))
    enddo ; enddo ; endif
    if (CS%id_KhTr_v > 0) then ; do J=js-1,je ;  do i=is,ie
      Kh_v(i,J) = CS%KhTr
      khdt_y(i,J) = dt*(CS%KhTr*(G%dx_v(i,J)*G%IDYv(i,J)))
    enddo ; enddo ; else ; do J=js-1,je ;  do i=is,ie
      khdt_y(i,J) = dt*(CS%KhTr*(G%dx_v(i,J)*G%IDYv(i,J)))
    enddo ; enddo ; endif
  endif

  if (CS%check_diffusive_CFL) then
    max_CFL = 0.0
    do j=js,je ; do i=is,ie
      CFL(i,j) = 2.0*((khdt_x(I-1,j) + khdt_x(I,j)) + &
                      (khdt_y(i,J-1) + khdt_y(i,J))) * G%IDXDYh(i,j)
      if (max_CFL < CFL(i,j)) max_CFL = CFL(i,j)
    enddo ; enddo
    call cpu_clock_begin(id_clock_sync)
    call max_across_PEs(max_CFL)
    call cpu_clock_end(id_clock_sync)
    num_itts = max(1,ceiling(max_CFL))
    I_numitts = 1.0 ; if (num_itts > 1) I_numitts = 1.0 / (real(num_itts))
  else
    num_itts = 1 ; I_numitts = 1.0
  endif

  do m=1,ntr
    if (associated(Tr(m)%df_x)) then
      do k=1,nz ; do j=js,je ; do I=is-1,ie
        Tr(m)%df_x(I,j,k) = 0.0
      enddo ; enddo ; enddo
    endif
    if (associated(Tr(m)%df_y)) then
      do k=1,nz ; do J=js-1,je ; do i=is,ie
        Tr(m)%df_y(i,J,k) = 0.0
      enddo ; enddo ; enddo
    endif
    if (associated(Tr(m)%df2d_x)) then
      do j=js,je ; do I=is-1,ie ; Tr(m)%df2d_x(I,j) = 0.0 ; enddo ; enddo
    endif
    if (associated(Tr(m)%df2d_y)) then
      do J=js-1,je ; do i=is,ie ; Tr(m)%df2d_y(i,J) = 0.0 ; enddo ; enddo
    endif
  enddo

  do itt=1,num_itts
    call cpu_clock_begin(id_clock_pass)
    do m=1,ntr-1
      call pass_var(Tr(m)%t(:,:,:), G%Domain, complete=.false.)
    enddo
    call pass_var(Tr(ntr)%t(:,:,:), G%Domain)
    call cpu_clock_end(id_clock_pass)

    do k=1,nz
      scale = I_numitts
      if (CS%Diffuse_ML_interior) then
        if (k<=G%nkml) then
          if (CS%ML_KhTr_scale <= 0.0) cycle
          scale = I_numitts * CS%ML_KhTr_scale
        endif
        if ((k>G%nkml) .and. (k<=G%nk_rho_varies)) cycle
      endif

      do J=js-1,je ; do i=is,ie
        Coef_y(i,J) = ((scale * khdt_y(i,J))*2.0*(h(i,j,k)*h(i,j+1,k))) / &
                                                 (h(i,j,k)+h(i,j+1,k)+h_neglect)
      enddo ; enddo

      do j=js,je
        do I=is-1,ie
          Coef_x(I,j) = ((scale * khdt_x(I,j))*2.0*(h(i,j,k)*h(i+1,j,k))) / &
                                                   (h(i,j,k)+h(i+1,j,k)+h_neglect)
        enddo

        do i=is,ie
          Ihdxdy(i,j) = G%IDXDYh(i,j) / (h(i,j,k)+h_neglect)
        enddo
      enddo

      do m=1,ntr
        do j=js,je ; do i=is,ie
          dTr(i,j) = Ihdxdy(i,j) * &
            ((Coef_x(I-1,j) * (Tr(m)%t(i-1,j,k) - Tr(m)%t(i,j,k)) - &
              Coef_x(I,j) * (Tr(m)%t(i,j,k) - Tr(m)%t(i+1,j,k))) + &
             (Coef_y(i,J-1) * (Tr(m)%t(i,j-1,k) - Tr(m)%t(i,j,k)) - &
              Coef_y(i,J) * (Tr(m)%t(i,j,k) - Tr(m)%t(i,j+1,k))))
        enddo ; enddo
        if (associated(Tr(m)%df_x)) then ; do j=js,je ; do I=G%Iscq,G%Iecq
          Tr(m)%df_x(I,j,k) = Tr(m)%df_x(I,j,k) + Coef_x(I,j) * &
                              (Tr(m)%t(i,j,k) - Tr(m)%t(i+1,j,k))*Idt
        enddo ; enddo ; endif
        if (associated(Tr(m)%df_y)) then ; do J=G%Jscq,G%Jecq ; do i=is,ie
          Tr(m)%df_y(i,J,k) = Tr(m)%df_y(i,J,k) + Coef_y(i,J) * &
                              (Tr(m)%t(i,j,k) - Tr(m)%t(i,j+1,k))*Idt
        enddo ; enddo ; endif
        if (associated(Tr(m)%df2d_x)) then ; do j=js,je ; do I=G%Iscq,G%Iecq
          Tr(m)%df2d_x(I,j) = Tr(m)%df2d_x(I,j) + Coef_x(I,j) * &
                              (Tr(m)%t(i,j,k) - Tr(m)%t(i+1,j,k))*Idt
        enddo ; enddo ; endif
        if (associated(Tr(m)%df2d_y)) then ; do J=G%Jscq,G%Jecq ; do i=is,ie
          Tr(m)%df2d_y(i,J) = Tr(m)%df2d_y(i,J) + Coef_y(i,J) * &
                              (Tr(m)%t(i,j,k) - Tr(m)%t(i,j+1,k))*Idt
        enddo ; enddo ; endif
        do j=js,je ; do i=is,ie
          Tr(m)%t(i,j,k) = Tr(m)%t(i,j,k) + dTr(i,j)
        enddo ; enddo
      enddo

    enddo ! End of k loop.

  enddo ! End of "while" loop.
  call cpu_clock_end(id_clock_diffuse)

  if (CS%Diffuse_ML_interior) then
    if (CS%debug) call MOM_tracer_chksum("Before epipycnal diff ", Tr, ntr, G)

    call cpu_clock_begin(id_clock_epimix)
    call tracer_epipycnal_ML_diff(h, dt, Tr, khdt_x, khdt_y, G, CS, tv, num_itts)
    call cpu_clock_end(id_clock_epimix)
  endif

  if (CS%debug) call MOM_tracer_chksum("After tracer diffusion ", Tr, ntr, G)

  if (CS%id_KhTr_u > 0) then
    do j=js,je ; do I=is-1,ie
      Kh_u(I,j) = G%umask(I,j)*Kh_u(I,j)
    enddo ; enddo
    call post_data(CS%id_KhTr_u, Kh_u, CS%diag)
  endif
  if (CS%id_KhTr_v > 0) then
    do J=js-1,je ; do i=is,ie
      Kh_v(i,J) = G%vmask(i,J)*Kh_v(i,J)
    enddo ; enddo
    call post_data(CS%id_KhTr_v, Kh_v, CS%diag)
  endif

end subroutine tracer_hordiff

subroutine tracer_epipycnal_ML_diff(h, dt, Tr, khdt_epi_x, khdt_epi_y, G, CS, &
                                    tv, num_itts)
  real, dimension(NXMEM_,NYMEM_,NKMEM_),  intent(in)    :: h
  real,                                intent(in)    :: dt
  type(tracer),                        intent(inout) :: Tr(:)
  real, dimension(NXMEMQ_,NYMEM_),     intent(in)    :: khdt_epi_x
  real, dimension(NXMEM_,NYMEMQ_),     intent(in)    :: khdt_epi_y
  type(ocean_grid_type),               intent(inout) :: G
  type(advect_tracer_CS),              intent(in)    :: CS
  type(thermo_var_ptrs),               intent(in)    :: tv
  integer,                             intent(in)    :: num_itts
!   This subroutine does epipycnal diffusion of all tracers between the mixed
! and buffer layers and the interior, using the diffusivity in CS%KhTr.
! Multiple iterations are used (if necessary) so that there is no limit on the
! acceptable time increment.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (in)      dt - Time increment in s.
!  (in/out)  Tr - An array of all of the registered tracers.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 advect_tracer_init.
!  (in/out?) tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs,
!                 and these may (probably will) point to some of the same arrays
!                 as Tr does.  tv is required for epipycnal mixing between the
!                 mixed layer and the interior.
!  (in)      num_itts - The number of iterations to use, usually 1.
  
  real, dimension(SZI_(G), SZJ_(G)) :: &
    Rml_max  ! The maximum coordinate density within the mixed layer, in kg m-3.
  real, dimension(SZI_(G), SZJ_(G), max(1,G%nk_rho_varies)) :: &
    rho_coord ! The coordinate density that is used to mix along, in kg m-3.

  ! The naming mnemnonic is a=above,b=below,L=Left,R=Right,u=u-point,v=v-point.
  ! These are 1-D arrays of pointers to 2-d arrays to minimize memory usage.
  type(p2d), dimension(SZJ_(G)) :: &
    deep_wt_Lu, deep_wt_Ru, &  ! The relative weighting of the deeper of a pair, ND.
    hP_Lu, hP_Ru       ! The total thickness on each side for each pair, in m or kg m-2.
    
  type(p2d), dimension(SZJQ_(G)) :: &
    deep_wt_Lv, deep_wt_Rv, & ! The relative weighting of the deeper of a pair, ND.
    hP_Lv, hP_Rv       ! The total thickness on each side for each pair, in m or kg m-2.

  type(p2di), dimension(SZJ_(G)) :: &
    k0b_Lu, k0a_Lu, &  ! The original k-indices of the layers that participate
    k0b_Ru, k0a_Ru     ! in each pair of mixing at u-faces.
  type(p2di), dimension(SZJQ_(G)) :: &
    k0b_Lv, k0a_Lv, &  ! The original k-indices of the layers that participate
    k0b_Rv, k0a_Rv     ! in each pair of mixing at v-faces.

  real, dimension(SZI_(G), SZJ_(G), SZK_(G)) :: &
    tr_flux_conv  ! The flux convergence of tracers, in TR m3 or TR kg.
 
  real, dimension(SZI_(G), SZK_(G), SZJ_(G)) :: &
    rho_srt, & ! The density of each layer of the sorted columns, in kg m-3.
    h_srt      ! The thickness of each layer of the sorted columns, in m or kg m-2.
  integer, dimension(SZI_(G), SZK_(G), SZJ_(G)) :: &
    k0_srt     ! The original k-index that each layer of the sorted column
               ! corresponds to.

  real, dimension(SZK_(G)) :: &
    h_demand_L, & ! The thickness in the left (_L) or right (_R) column that
    h_demand_R, & ! is demanded to match the thickness in the counterpart, in H.
    h_used_L, &   ! The summed thickness from the left or right columns that
    h_used_R, &   ! have actually been used, in m or kg m-2 (H).
    h_supply_frac_L, &  ! The fraction of the demanded thickness that can
    h_supply_frac_R     ! actually be supplied from a layer.
  integer, dimension(SZK_(G)) :: &
    kbs_Lp, &   ! The sorted indicies of the Left and Right columns for
    kbs_Rp      ! each pairing.

  integer, dimension(SZI_(G), SZJ_(G))  :: &
    num_srt, &   ! The number of layers that are sorted in each column.
    k_end_srt, & ! The maximum index in each column that might need to be
                 ! sorted, based on neighboring values of max_kRho
    max_kRho     ! The index of the layer whose target density is just denser
                 ! than the densest part of the mixed layer.
  integer, dimension(SZJ_(G))           :: &
    max_srt      ! The maximum value of num_srt in a k-row.
  integer, dimension(SZIQ_(G), SZJ_(G)) :: &
    nPu          ! The number of epipycnal pairings at each u-point.
  integer, dimension(SZI_(G), SZJQ_(G)) :: &
    nPv          ! The number of epipycnal pairings at each v-point.
  real :: h_exclude    ! A thickness that layers must attain to be considered
                       ! for inclusion in mixing, in m.
  real :: Idt        ! The inverse of the time step, in s-1.
  real :: I_maxitt   ! The inverse of the maximum number of iterations.
  real :: rho_pair, rho_a, rho_b  ! Temporary densities, in kg m-3.
  real :: Tr_min_face  ! The minimum and maximum tracer concentrations
  real :: Tr_max_face  ! associated with a pairing, in conc.
  real :: Tr_La, Tr_Lb ! The 4 tracer concentrations that might be
  real :: Tr_Ra, Tr_Rb ! associated with a pairing, in conc.
  real :: Tr_av_L    ! The average tracer concentrations on the left and right
  real :: Tr_av_R    ! sides of a pairing, in conc.
  real :: Tr_flux    ! The tracer flux from left to right in a pair, in conc m3.
  real :: Tr_adj_vert  ! A downward vertical adjustment to Tr_flux between the
                     ! two cells that make up one side of the pairing, in conc m3.
  real :: h_L, h_R   ! Thicknesses to the left and right, in m or kg m-2 (H).
  real :: wt_a, wt_b ! Fractional weights of layers above and below, ND.
  real :: vol        ! A cell volume or mass, in m3 or kg (H m2).
  logical, dimension(SZK_(G)) :: &
    left_set, &  ! If true, the left or right point determines the density of
    right_set    ! of the trio.  If densities are exactly equal, both are true.
  real :: tmp
  real :: p_ref_cv(SZI_(G))

  integer :: k_max, k_min, k_test, itmp
  integer :: i, j, k, k2, m, is, ie, js, je, nz, ntr, nkmb
  integer :: isd, ied, jsd, jed, Isdq, Iedq, k_size
  integer :: kL, kR, kLa, kLb, kRa, kRb, nP, itt, ns, max_itt
  integer :: PEmax_kRho
  integer :: isv, iev, jsv, jev ! The valid range of the indices.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isdq = G%Isdq ; Iedq = G%Iedq
  ntr = CS%ntr
  Idt = 1.0/dt
  nkmb = G%nk_rho_varies

  if (num_itts <= 1) then
    max_itt = 1 ; I_maxitt = 1.0
  else
    max_itt = num_itts ; I_maxitt = 1.0 / (real(max_itt))
  endif

  do i=is-2,ie+2 ; p_ref_cv(i) = tv%P_Ref ; enddo

  call cpu_clock_begin(id_clock_pass)
  do m=1,ntr-1
    call pass_var(Tr(m)%t(:,:,:), G%Domain, complete=.false.)
  enddo
  call pass_var(Tr(ntr)%t(:,:,:), G%Domain)
  call cpu_clock_end(id_clock_pass)

  ! Determine which layers the mixed- and buffer-layers map into...
  do k=1,nkmb ; do j=js-2,je+2
    call calculate_density(tv%T(:,j,k),tv%S(:,j,k), p_ref_cv, &
                           rho_coord(:,j,k), is-2, ie-is+5, tv%eqn_of_state)
  enddo ; enddo

  do j=js-2,je+2 ; do i=is-2,ie+2
    Rml_max(i,j) = rho_coord(i,j,1)
    num_srt(i,j) = 0 ; max_kRho(i,j) = 0
  enddo ; enddo
  do k=2,nkmb ; do j=js-2,je+2 ; do i=is-2,ie+2
    if (Rml_max(i,j) < rho_coord(i,j,k)) Rml_max(i,j) = rho_coord(i,j,k)
  enddo ; enddo ; enddo

  !   Use bracketing and bisection to find the k-level that the densest of the
  ! mixed and buffer layer corresponds to, such that:
  !     G%Rlay(max_kRho-1) < Rml_max <= G%Rlay(max_kRho)
  do j=js-2,je+2 ; do i=is-2,ie+2 ; if (G%hmask(i,j) > 0.5) then
    if (Rml_max(i,j) > G%Rlay(nz)) then ; max_kRho(i,j) = nz+1
    elseif (Rml_max(i,j) <= G%Rlay(nkmb+1)) then ; max_kRho(i,j) = nkmb+1
    else
      k_min = nkmb+2 ; k_max = nz
      do
        k_test = (k_min + k_max) / 2
        if (Rml_max(i,j) <= G%Rlay(k_test-1)) then ; k_max = k_test-1
        elseif (G%Rlay(k_test) < Rml_max(i,j)) then ; k_min = k_test+1
        else ; max_kRho(i,j) = k_test ; exit ; endif

        if (k_min == k_max) then ; max_kRho(i,j) = k_max ; exit ; endif
      enddo
    endif
  endif ; enddo ; enddo

  PEmax_kRho = 0
  do j=js-1,je+1 ; do i=is-1,ie+1
    k_end_srt(i,j) = max(max_kRho(i,j), max_kRho(i-1,j), max_kRho(i+1,j), &
                         max_kRho(i,j-1), max_kRho(i,j+1))
    if (PEmax_kRho < k_end_srt(i,j)) PEmax_kRho = k_end_srt(i,j)
  enddo ; enddo
  if (PEmax_kRho > nz) PEmax_kRho = nz ! PEmax_kRho could have been nz+1.

  h_exclude = 10.0*(G%Angstrom + G%H_subroundoff)

  do k=1,nkmb ; do j=js-1,je+1 ; do i=is-1,ie+1 ; if (G%hmask(i,j) > 0.5) then
    if (h(i,j,k) > h_exclude) then
      num_srt(i,j) = num_srt(i,j) + 1 ; ns = num_srt(i,j)
      k0_srt(i,ns,j) = k
      rho_srt(i,ns,j) = rho_coord(i,j,k)
      h_srt(i,ns,j) = h(i,j,k)
    endif
  endif ; enddo ; enddo ; enddo
  do k=nkmb+1,PEmax_kRho ; do j=js-1,je+1 ; do i=is-1,ie+1 ; if (G%hmask(i,j) > 0.5) then
    if ((k<=k_end_srt(i,j)) .and. (h(i,j,k) > h_exclude)) then
      num_srt(i,j) = num_srt(i,j) + 1 ; ns = num_srt(i,j)
      k0_srt(i,ns,j) = k
      rho_srt(i,ns,j) = G%Rlay(k)
      h_srt(i,ns,j) = h(i,j,k)
    endif
  endif ; enddo ; enddo ; enddo

  ! Sort each column by increasing density.  This should already be close, 
  ! and the size of the arrays are small, so straight insertion is used.
  do j=js-1,je+1 ; do i=is-1,ie+1
    do k=2,num_srt(i,j) ; if (rho_srt(i,k,j) < rho_srt(i,k-1,j)) then
      ! The last segment needs to be shuffled earlier in the list.
      do k2 = k,2,-1 ; if (rho_srt(i,k2,j) >= rho_srt(i,k2-1,j)) exit
        itmp = k0_srt(i,k2-1,j) ; k0_srt(i,k2-1,j) = k0_srt(i,k2,j) ; k0_srt(i,k2,j) = itmp
        tmp = rho_srt(i,k2-1,j) ; rho_srt(i,k2-1,j) = rho_srt(i,k2,j) ; rho_srt(i,k2,j) = tmp
        tmp = h_srt(i,k2-1,j) ; h_srt(i,k2-1,j) = h_srt(i,k2,j) ; h_srt(i,k2,j) = tmp
      enddo
    endif ; enddo
  enddo ; enddo

  do j=js-1,je+1
    max_srt(j) = 0
    do i=is-1,ie+1 ; max_srt(j) = max(max_srt(j), num_srt(i,j)) ; enddo
  enddo

  do j=js,je
    k_size = max(2*max_srt(j),1)
    allocate(deep_wt_Lu(j)%p(Isdq:Iedq,k_size))
    allocate(deep_wt_Ru(j)%p(Isdq:Iedq,k_size))
    allocate(hP_Lu(j)%p(Isdq:Iedq,k_size))
    allocate(hP_Ru(j)%p(Isdq:Iedq,k_size))
    allocate(k0a_Lu(j)%p(Isdq:Iedq,k_size))
    allocate(k0a_Ru(j)%p(Isdq:Iedq,k_size))
    allocate(k0b_Lu(j)%p(Isdq:Iedq,k_size))
    allocate(k0b_Ru(j)%p(Isdq:Iedq,k_size))
  enddo

  do j=js,je ; do I=is-1,ie ; if (G%umask(I,j) > 0.5) then
    ! Set up the pairings for fluxes through the zonal faces.

    do k=1,num_srt(i,j)   ; h_demand_L(k) = 0.0 ; h_used_L(k) = 0.0 ; enddo
    do k=1,num_srt(i+1,j) ; h_demand_R(k) = 0.0 ; h_used_R(k) = 0.0 ; enddo

    ! First merge the left and right lists into a single, sorted list.

    !   Discard any layers that are lighter than the lightest in the other
    ! column.  They can only participate in mixing as the lighter part of a
    ! pair of points.
    if (rho_srt(i,1,j) < rho_srt(i+1,1,j)) then
      kR = 1
      do kL=2,num_srt(i,j) ; if (rho_srt(i,kL,j) >= rho_srt(i+1,1,j)) exit ; enddo
    elseif (rho_srt(i+1,1,j) < rho_srt(i,1,j)) then
      kL = 1
      do kR=2,num_srt(i+1,j) ; if (rho_srt(i+1,kR,j) >= rho_srt(i,1,j)) exit ; enddo
    else
      kL = 1 ; kR = 1
    endif
    nP = 0
    do ! Loop to accumulate pairs of columns.
      if ((kL > num_srt(i,j)) .or. (kR > num_srt(i+1,j))) exit

      if (rho_srt(i,kL,j) > rho_srt(i+1,kR,j)) then
      ! The right point is lighter and defines the density for this trio.
        nP = nP+1 ; k = nP
        rho_pair = rho_srt(i+1,kR,j)

        k0b_Lu(j)%p(I,k) = k0_srt(i,kL,j) ; k0b_Ru(j)%p(I,k) = k0_srt(i+1,kR,j)
        k0a_Lu(j)%p(I,k) = k0_srt(i,kL-1,j) ; k0a_Ru(j)%p(I,k) = k0b_Ru(j)%p(I,k)
        kbs_Lp(k) = kL ; kbs_Rp(k) = kR

        rho_a = rho_srt(i,kL-1,j) ; rho_b = rho_srt(i,kL,j)
        wt_b = 1.0 ; if (abs(rho_a - rho_b) > abs(rho_pair - rho_a)) &
          wt_b = (rho_pair - rho_a) / (rho_b - rho_a)
        deep_wt_Lu(j)%p(I,k) = wt_b ; deep_wt_Ru(j)%p(I,k) = 1.0

        h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i+1,kR,j) * wt_b
        h_demand_L(kL-1) = h_demand_L(kL-1) + 0.5*h_srt(i+1,kR,j) * (1.0-wt_b)

        kR = kR+1 ; left_set(k) = .false. ; right_set(k) = .true.
      elseif (rho_srt(i,kL,j) < rho_srt(i+1,kR,j)) then
      ! The left point is lighter and defines the density for this trio.
        nP = nP+1 ; k = nP
        rho_pair = rho_srt(i,kL,j)
        k0b_Lu(j)%p(I,k) = k0_srt(i,kL,j) ; k0b_Ru(j)%p(I,k) = k0_srt(i+1,kR,j)
        k0a_Lu(j)%p(I,k) = k0b_Lu(j)%p(I,k) ; k0a_Ru(j)%p(I,k) = k0_srt(i+1,kR-1,j)

        kbs_Lp(k) = kL ; kbs_Rp(k) = kR

        rho_a = rho_srt(i+1,kR-1,j) ; rho_b = rho_srt(i+1,kR,j)
        wt_b = 1.0 ; if (abs(rho_a - rho_b) > abs(rho_pair - rho_a)) &
          wt_b = (rho_pair - rho_a) / (rho_b - rho_a)
        deep_wt_Lu(j)%p(I,k) = 1.0 ; deep_wt_Ru(j)%p(I,k) = wt_b

        h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j) * wt_b
        h_demand_R(kR-1) = h_demand_R(kR-1) + 0.5*h_srt(i,kL,j) * (1.0-wt_b)

        kL = kL+1 ; left_set(k) = .true. ; right_set(k) = .false.
      elseif ((k0_srt(i,kL,j) <= nkmb) .or. (k0_srt(i+1,kR,j) <= nkmb)) then
        ! The densities are exactly equal and one layer is above the interior.
        nP = nP+1 ; k = nP
        k0b_Lu(j)%p(I,k) = k0_srt(i,kL,j) ; k0b_Ru(j)%p(I,k) = k0_srt(i+1,kR,j)
        k0a_Lu(j)%p(I,k) = k0b_Lu(j)%p(I,k) ; k0a_Ru(j)%p(I,k) = k0b_Ru(j)%p(I,k)
        kbs_Lp(k) = kL ; kbs_Rp(k) = kR
        deep_wt_Lu(j)%p(I,k) = 1.0 ; deep_wt_Ru(j)%p(I,k) = 1.0

        h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i+1,kR,j)
        h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j)

        kL = kL+1 ; kR = kR+1 ; left_set(k) = .true. ; right_set(k) = .true.
      else ! The densities are exactly equal and in the interior.
        ! Mixing in this case has already occurred, so accumulate the thickness
        ! demanded for that mixing and skip onward.
        h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i+1,kR,j) 
        h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j)

        kL = kL+1 ; kR = kR+1
      endif
    enddo ! Loop to accumulate pairs of columns.
    nPu(I,j) = nP ! This is the number of active pairings.

    ! Determine what fraction of the thickness "demand" can be supplied.
    do k=1,num_srt(i+1,j)
      h_supply_frac_R(k) = 1.0
      if (h_demand_R(k) > 0.5*h_srt(i+1,k,j)) &
        h_supply_frac_R(k) = 0.5*h_srt(i+1,k,j) / h_demand_R(k)
    enddo
    do k=1,num_srt(i,j)
      h_supply_frac_L(k) = 1.0
      if (h_demand_L(k) > 0.5*h_srt(i,k,j)) &
        h_supply_frac_L(k) = 0.5*h_srt(i,k,j) / h_demand_L(k)
    enddo

    !  Distribute the "exported" thicknesses proportionately.
    do k=1,nPu(I,j)
      kL = kbs_Lp(k) ; kR = kbs_Rp(k)
      hP_Lu(j)%p(I,k) = 0.0 ; hP_Ru(j)%p(I,k) = 0.0
      if (left_set(k)) then ! Add the contributing thicknesses on the right.
        if (deep_wt_Ru(j)%p(I,k) < 1.0) then
          hP_Ru(j)%p(I,k) = 0.5*h_srt(i,kL,j) * min(h_supply_frac_R(kR), h_supply_frac_R(kR-1))
          wt_b = deep_wt_Ru(j)%p(I,k)
          h_used_R(kR-1) = h_used_R(kR-1) + (1.0 - wt_b)*hP_Ru(j)%p(I,k)
          h_used_R(kR) = h_used_R(kR) + wt_b*hP_Ru(j)%p(I,k)
        else
          hP_Ru(j)%p(I,k) = 0.5*h_srt(i,kL,j) * h_supply_frac_R(kR)
          h_used_R(kR) = h_used_R(kR) + hP_Ru(j)%p(I,k)
        endif
      endif
      if (right_set(k)) then ! Add the contributing thicknesses on the left.
        if (deep_wt_Lu(j)%p(I,k) < 1.0) then
          hP_Lu(j)%p(I,k) = 0.5*h_srt(i+1,kR,j) * min(h_supply_frac_L(kL), h_supply_frac_L(kL-1))
          wt_b = deep_wt_Lu(j)%p(I,k)
          h_used_L(kL-1) = h_used_L(kL-1) + (1.0 - wt_b)*hP_Lu(j)%p(I,k)
          h_used_L(kL) = h_used_L(kL) + wt_b*hP_Lu(j)%p(I,k)
        else
          hP_Lu(j)%p(I,k) = 0.5*h_srt(i+1,kR,j) * h_supply_frac_L(kL)
          h_used_L(kL) = h_used_L(kL) + hP_Lu(j)%p(I,k)
        endif
      endif
    enddo

    !   The left-over thickness (at least half the layer thickness) is now
    ! added to the thicknesses of the importing columns.
    do k=1,nPu(I,j)
      if (left_set(k)) hP_Lu(j)%p(I,k) = hP_Lu(j)%p(I,k) + &
                           (h_srt(i,kbs_Lp(k),j) - h_used_L(kbs_Lp(k)))
      if (right_set(k)) hP_Ru(j)%p(I,k) = hP_Ru(j)%p(I,k) + &
                            (h_srt(i+1,kbs_Rp(k),j) - h_used_R(kbs_Rp(k)))
    enddo

  endif ; enddo ; enddo ! i- & j- loops over zonal faces.

  do J=js-1,je
    k_size = max(max_srt(j)+max_srt(j+1),1)
    allocate(deep_wt_Lv(J)%p(isd:ied,k_size))
    allocate(deep_wt_Rv(J)%p(isd:ied,k_size))
    allocate(hP_Lv(J)%p(isd:ied,k_size))
    allocate(hP_Rv(J)%p(isd:ied,k_size))
    allocate(k0a_Lv(J)%p(isd:ied,k_size))
    allocate(k0a_Rv(J)%p(isd:ied,k_size))
    allocate(k0b_Lv(J)%p(isd:ied,k_size))
    allocate(k0b_Rv(J)%p(isd:ied,k_size))
  enddo

  do J=js-1,je ; do i=is,ie ; if (G%vmask(i,J) > 0.5) then
    ! Set up the pairings for fluxes through the meridional faces.

    do k=1,num_srt(i,j)   ; h_demand_L(k) = 0.0 ; h_used_L(k) = 0.0 ; enddo
    do k=1,num_srt(i,j+1) ; h_demand_R(k) = 0.0 ; h_used_R(k) = 0.0 ; enddo

    ! First merge the left and right lists into a single, sorted list.

    !   Discard any layers that are lighter than the lightest in the other
    ! column.  They can only participate in mixing as the lighter part of a
    ! pair of points.
    if (rho_srt(i,1,j) < rho_srt(i,1,j+1)) then
      kR = 1
      do kL=2,num_srt(i,j) ; if (rho_srt(i,kL,j) >= rho_srt(i,1,j+1)) exit ; enddo
    elseif (rho_srt(i,1,j+1) < rho_srt(i,1,j)) then
      kL = 1
      do kR=2,num_srt(i,j+1) ; if (rho_srt(i,kR,j+1) >= rho_srt(i,1,j)) exit ; enddo
    else
      kL = 1 ; kR = 1
    endif
    nP = 0
    do ! Loop to accumulate pairs of columns.
      if ((kL > num_srt(i,j)) .or. (kR > num_srt(i,j+1))) exit

      if (rho_srt(i,kL,j) > rho_srt(i,kR,j+1)) then
      ! The right point is lighter and defines the density for this trio.
        nP = nP+1 ; k = nP
        rho_pair = rho_srt(i,kR,j+1)

        k0b_Lv(J)%p(i,k) = k0_srt(i,kL,j)   ; k0b_Rv(J)%p(i,k) = k0_srt(i,kR,j+1)
        k0a_Lv(J)%p(i,k) = k0_srt(i,kL-1,j) ; k0a_Rv(J)%p(i,k) = k0b_Rv(J)%p(i,k)
        kbs_Lp(k) = kL ; kbs_Rp(k) = kR

        rho_a = rho_srt(i,kL-1,j) ; rho_b = rho_srt(i,kL,j)
        wt_b = 1.0 ; if (abs(rho_a - rho_b) > abs(rho_pair - rho_a)) &
          wt_b = (rho_pair - rho_a) / (rho_b - rho_a)
        deep_wt_Lv(J)%p(i,k) = wt_b ; deep_wt_Rv(J)%p(i,k) = 1.0

        h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i,kR,j+1) * wt_b
        h_demand_L(kL-1) = h_demand_L(kL-1) + 0.5*h_srt(i,kR,j+1) * (1.0-wt_b)

        kR = kR+1 ; left_set(k) = .false. ; right_set(k) = .true.
      elseif (rho_srt(i,kL,j) < rho_srt(i,kR,j+1)) then
      ! The left point is lighter and defines the density for this trio.
        nP = nP+1 ; k = nP
        rho_pair = rho_srt(i,kL,j)
        k0b_Lv(J)%p(i,k) = k0_srt(i,kL,j) ; k0b_Rv(J)%p(i,k) = k0_srt(i,kR,j+1)
        k0a_Lv(J)%p(i,k) = k0b_Lv(J)%p(i,k) ; k0a_Rv(J)%p(i,k) = k0_srt(i,kR-1,j+1)

        kbs_Lp(k) = kL ; kbs_Rp(k) = kR

        rho_a = rho_srt(i,kR-1,j+1) ; rho_b = rho_srt(i,kR,j+1)
        wt_b = 1.0 ; if (abs(rho_a - rho_b) > abs(rho_pair - rho_a)) &
          wt_b = (rho_pair - rho_a) / (rho_b - rho_a)
        deep_wt_Lv(J)%p(i,k) = 1.0 ; deep_wt_Rv(J)%p(i,k) = wt_b

        h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j) * wt_b
        h_demand_R(kR-1) = h_demand_R(kR-1) + 0.5*h_srt(i,kL,j) * (1.0-wt_b)

        kL = kL+1 ; left_set(k) = .true. ; right_set(k) = .false.
      elseif ((k0_srt(i,kL,j) <= nkmb) .or. (k0_srt(i,kR,j+1) <= nkmb)) then
        ! The densities are exactly equal and one layer is above the interior.
        nP = nP+1 ; k = nP
        k0b_Lv(J)%p(i,k) = k0_srt(i,kL,j) ; k0b_Rv(J)%p(i,k) = k0_srt(i,kR,j+1)
        k0a_Lv(J)%p(i,k) = k0b_Lv(J)%p(i,k)  ; k0a_Rv(J)%p(i,k) = k0b_Rv(J)%p(i,k)
        kbs_Lp(k) = kL ; kbs_Rp(k) = kR
        deep_wt_Lv(J)%p(i,k) = 1.0 ; deep_wt_Rv(J)%p(i,k) = 1.0

        h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i,kR,j+1)
        h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j)

        kL = kL+1 ; kR = kR+1 ; left_set(k) = .true. ; right_set(k) = .true.
      else ! The densities are exactly equal and in the interior.
        ! Mixing in this case has already occurred, so accumulate the thickness
        ! demanded for that mixing and skip onward.
        h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i,kR,j+1) 
        h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j)

        kL = kL+1 ; kR = kR+1
      endif
    enddo ! Loop to accumulate pairs of columns.
    nPv(i,J) = nP ! This is the number of active pairings.

    ! Determine what fraction of the thickness "demand" can be supplied.
    do k=1,num_srt(i,j+1)
      h_supply_frac_R(k) = 1.0
      if (h_demand_R(k) > 0.5*h_srt(i,k,j+1)) &
        h_supply_frac_R(k) = 0.5*h_srt(i,k,j+1) / h_demand_R(k)
    enddo
    do k=1,num_srt(i,j)
      h_supply_frac_L(k) = 1.0
      if (h_demand_L(k) > 0.5*h_srt(i,k,j)) &
        h_supply_frac_L(k) = 0.5*h_srt(i,k,j) / h_demand_L(k)
    enddo

    !  Distribute the "exported" thicknesses proportionately.
    do k=1,nPv(i,J)
      kL = kbs_Lp(k) ; kR = kbs_Rp(k)
      hP_Lv(J)%p(i,k) = 0.0 ; hP_Rv(J)%p(i,k) = 0.0
      if (left_set(k)) then ! Add the contributing thicknesses on the right.
        if (deep_wt_Rv(J)%p(i,k) < 1.0) then
          hP_Rv(J)%p(i,k) = 0.5*h_srt(i,kL,j) * min(h_supply_frac_R(kR), h_supply_frac_R(kR-1))
          wt_b = deep_wt_Rv(J)%p(i,k)
          h_used_R(kR-1) = h_used_R(kR-1) + (1.0 - wt_b) * hP_Rv(J)%p(i,k)
          h_used_R(kR) = h_used_R(kR) + wt_b * hP_Rv(J)%p(i,k)
        else
          hP_Rv(J)%p(i,k) = 0.5*h_srt(i,kL,j) * h_supply_frac_R(kR)
          h_used_R(kR) = h_used_R(kR) + hP_Rv(J)%p(i,k)
        endif
      endif
      if (right_set(k)) then ! Add the contributing thicknesses on the left.
        if (deep_wt_Lv(J)%p(i,k) < 1.0) then
          hP_Lv(J)%p(i,k) = 0.5*h_srt(i,kR,j+1) * min(h_supply_frac_L(kL), h_supply_frac_L(kL-1))
          wt_b = deep_wt_Lv(J)%p(i,k)
          h_used_L(kL-1) = h_used_L(kL-1) + (1.0 - wt_b) * hP_Lv(J)%p(i,k)
          h_used_L(kL) = h_used_L(kL) + wt_b * hP_Lv(J)%p(i,k)
        else
          hP_Lv(J)%p(i,k) = 0.5*h_srt(i,kR,j+1) * h_supply_frac_L(kL)
          h_used_L(kL) = h_used_L(kL) + hP_Lv(J)%p(i,k)
        endif
      endif
    enddo

    !   The left-over thickness (at least half the layer thickness) is now
    ! added to the thicknesses of the importing columns.
    do k=1,nPv(i,J)
      if (left_set(k)) hP_Lv(J)%p(i,k) = hP_Lv(J)%p(i,k) + &
                            (h_srt(i,kbs_Lp(k),j) - h_used_L(kbs_Lp(k)))
      if (right_set(k)) hP_Rv(J)%p(i,k) = hP_Rv(J)%p(i,k) + &
                             (h_srt(i,kbs_Rp(k),j+1) - h_used_R(kbs_Rp(k)))
    enddo


  endif ; enddo ; enddo ! i- & j- loops over meridional faces.

! The tracer-specific calculations start here.

  ! Zero out tracer tendencies.
  do k=1,PEmax_kRho ; do j=js-1,je+1 ; do i=is-1,ie+1
    tr_flux_conv(i,j,k) = 0.0
  enddo ; enddo ; enddo

  do itt=1,max_itt

    if (itt > 1) then ! The halos have already been filled if itt==1.
      call cpu_clock_begin(id_clock_pass)
      do m=1,ntr-1
        call pass_var(Tr(m)%t(:,:,:), G%Domain, complete=.false.)
      enddo
      call pass_var(Tr(ntr)%t(:,:,:), G%Domain)
      call cpu_clock_end(id_clock_pass)
    endif

    do m=1,ntr
    
      do j=js,je ; do I=is-1,ie ; if (G%umask(I,j) > 0.5) then
        ! Determine the fluxes through the zonal faces.

        ! Find the acceptable range of tracer concentration around this face.
        if (nPu(I,j) >= 1) then
          Tr_min_face = min(Tr(m)%t(i,j,1), Tr(m)%t(i+1,j,1))
          Tr_max_face = max(Tr(m)%t(i,j,1), Tr(m)%t(i+1,j,1))
          do k=2,nkmb
            Tr_min_face = min(Tr_min_face, Tr(m)%t(i,j,k), Tr(m)%t(i+1,j,k))
            Tr_max_face = max(Tr_max_face, Tr(m)%t(i,j,k), Tr(m)%t(i+1,j,k))
          enddo

          ! Include the next two layers denser than the densest buffer layer.
          kLa = nkmb+1 ; if (max_kRho(i,j) < nz+1) kLa = max_kRho(i,j)
          kLb = kLa ; if (max_kRho(i,j) < nz) kLb = max_kRho(i,j)+1
          kRa = nkmb+1 ; if (max_kRho(i+1,j) < nz+1) kRa = max_kRho(i+1,j)
          kRb = kRa ; if (max_kRho(i+1,j) < nz) kRb = max_kRho(i+1,j)+1
          Tr_La = Tr_min_face ; Tr_Lb = Tr_La ; Tr_Ra = Tr_La ; Tr_Rb = Tr_La
          if (h(i,j,kLa) > h_exclude) Tr_La = Tr(m)%t(i,j,kLa)
          if (h(i,j,kLb) > h_exclude) Tr_La = Tr(m)%t(i,j,kLb)
          if (h(i+1,j,kRa) > h_exclude) Tr_Ra = Tr(m)%t(i+1,j,kRa)
          if (h(i+1,j,kRb) > h_exclude) Tr_Rb = Tr(m)%t(i+1,j,kRb)
          Tr_min_face = min(Tr_min_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
          Tr_max_face = max(Tr_max_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)

          ! Include all points in diffusive pairings at this face.
          do k=1,nPu(I,j)
            Tr_Lb = Tr(m)%t(i,j,k0b_Lu(j)%p(I,k))
            Tr_Rb = Tr(m)%t(i+1,j,k0b_Ru(j)%p(I,k))
            Tr_La = Tr_Lb ; Tr_Ra = Tr_Rb
            if (deep_wt_Lu(j)%p(I,k) < 1.0) Tr_La = Tr(m)%t(i,j,k0a_Lu(j)%p(I,k))
            if (deep_wt_Ru(j)%p(I,k) < 1.0) Tr_Ra = Tr(m)%t(i+1,j,k0a_Ru(j)%p(I,k))
            Tr_min_face = min(Tr_min_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
            Tr_max_face = max(Tr_max_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
          enddo
        endif

        do k=1,nPu(I,j)
          kLb = k0b_Lu(j)%p(I,k) ; Tr_Lb = Tr(m)%t(i,j,kLb) ; Tr_av_L = Tr_Lb
          if (deep_wt_Lu(j)%p(I,k) < 1.0) then
            kLa = k0a_Lu(j)%p(I,k) ; Tr_La = Tr(m)%t(i,j,kLa)
            wt_b = deep_wt_Lu(j)%p(I,k)
            Tr_av_L = wt_b*Tr_Lb + (1.0-wt_b)*Tr_La
          endif

          kRb = k0b_Ru(j)%p(I,k) ; Tr_Rb = Tr(m)%t(i+1,j,kRb) ; Tr_av_R = Tr_Rb
          if (deep_wt_Ru(j)%p(I,k) < 1.0) then
            kRa = k0a_Ru(j)%p(I,k) ; Tr_Ra = Tr(m)%t(i+1,j,kRa)
            wt_b = deep_wt_Ru(j)%p(I,k)
            Tr_av_R = wt_b*Tr_Rb + (1.0-wt_b)*Tr_Ra
          endif

          h_L = hP_Lu(j)%p(I,k) ; h_R = hP_Ru(j)%p(I,k)
          Tr_flux = I_maxitt * khdt_epi_x(I,j) * (Tr_av_L - Tr_av_R) * &
            ((2.0 * h_L * h_R) / (h_L + h_R))
                    

          if (deep_wt_Lu(j)%p(I,k) >= 1.0) then
            tr_flux_conv(i,j,kLb) = tr_flux_conv(i,j,kLb) - Tr_flux
          else
            Tr_adj_vert = 0.0
            wt_b = deep_wt_Lu(j)%p(I,k) ; wt_a = 1.0 - wt_b
            vol = hP_Lu(j)%p(I,k) * G%DXDYh(i,j)

            !   Ensure that the tracer flux does not drive the tracer values
            ! outside of the range Tr_min_face <= Tr <= Tr_max_face, or if it
            ! does that the concentration in both contributing peices exceed
            ! this range equally. With downgradient fluxes and the initial tracer
            ! concentrations determining the valid range, the latter condition
            ! only enters for large values of the effective diffusive CFL number.
            if (Tr_flux > 0.0) then
              if (Tr_La < Tr_Lb) then ; if (vol*(Tr_La-Tr_min_face) < Tr_flux) &
                Tr_adj_vert = -wt_a * min(Tr_flux - vol * (Tr_La-Tr_min_face), &
                                          (vol*wt_b) * (Tr_Lb - Tr_La))
              else ; if (vol*(Tr_Lb-Tr_min_face) < Tr_flux) &
                Tr_adj_vert = wt_b * min(Tr_flux - vol * (Tr_Lb-Tr_min_face), &
                                         (vol*wt_a) * (Tr_La - Tr_Lb))
              endif
            elseif (Tr_flux < 0.0) then
              if (Tr_La > Tr_Lb) then ; if (vol * (Tr_max_face-Tr_La) < -Tr_flux) &
                Tr_adj_vert = wt_a * min(-Tr_flux - vol * (Tr_max_face-Tr_La), &
                                         (vol*wt_b) * (Tr_La - Tr_Lb))
              else ; if (vol*(Tr_max_face-Tr_Lb) < -Tr_flux) &
                Tr_adj_vert = -wt_b * min(-Tr_flux - vol * (Tr_max_face-Tr_Lb), &
                                          (vol*wt_a)*(Tr_Lb - Tr_La))
              endif
            endif

            tr_flux_conv(i,j,kLa) = tr_flux_conv(i,j,kLa) - (wt_a*Tr_flux + Tr_adj_vert)
            tr_flux_conv(i,j,kLb) = tr_flux_conv(i,j,kLb) - (wt_b*Tr_flux - Tr_adj_vert)
          endif

          if (deep_wt_Ru(j)%p(I,k) >= 1.0) then
            tr_flux_conv(i+1,j,kRb) = tr_flux_conv(i+1,j,kRb) + Tr_flux
          else
            Tr_adj_vert = 0.0
            wt_b = deep_wt_Ru(j)%p(I,k) ; wt_a = 1.0 - wt_b
            vol = hP_Ru(j)%p(I,k) * G%DXDYh(i+1,j)

            !   Ensure that the tracer flux does not drive the tracer values
            ! outside of the range Tr_min_face <= Tr <= Tr_max_face, or if it
            ! does that the concentration in both contributing peices exceed
            ! this range equally. With downgradient fluxes and the initial tracer
            ! concentrations determining the valid range, the latter condition
            ! only enters for large values of the effective diffusive CFL number.
            if (Tr_flux < 0.0) then
              if (Tr_Ra < Tr_Rb) then ; if (vol * (Tr_Ra-Tr_min_face) < -Tr_flux) &
                Tr_adj_vert = -wt_a * min(-Tr_flux - vol * (Tr_Ra-Tr_min_face), &
                                          (vol*wt_b) * (Tr_Rb - Tr_Ra))
              else ; if (vol*(Tr_Rb-Tr_min_face) < (-Tr_flux)) &
                Tr_adj_vert = wt_b * min(-Tr_flux - vol * (Tr_Rb-Tr_min_face), &
                                         (vol*wt_a) * (Tr_Ra - Tr_Rb))
              endif
            elseif (Tr_flux > 0.0) then
              if (Tr_Ra > Tr_Rb) then ; if (vol * (Tr_max_face-Tr_Ra) < Tr_flux) &
                Tr_adj_vert = wt_a * min(Tr_flux - vol * (Tr_max_face-Tr_Ra), &
                                         (vol*wt_b) * (Tr_Ra - Tr_Rb))
              else ; if (vol*(Tr_max_face-Tr_Rb) < Tr_flux) &
                Tr_adj_vert = -wt_b * min(Tr_flux - vol * (Tr_max_face-Tr_Rb), &
                                          (vol*wt_a)*(Tr_Rb - Tr_Ra))
              endif
            endif

            tr_flux_conv(i+1,j,kRa) = tr_flux_conv(i+1,j,kRa) + &
                                            (wt_a*Tr_flux - Tr_adj_vert)
            tr_flux_conv(i+1,j,kRb) = tr_flux_conv(i+1,j,kRb) + &
                                            (wt_b*Tr_flux + Tr_adj_vert)
          endif
          if (associated(Tr(m)%df2d_x)) &
            Tr(m)%df2d_x(I,j) = Tr(m)%df2d_x(I,j) + Tr_flux * Idt
        enddo ! Loop over pairings at faces.
      endif ; enddo ; enddo ! i- & j- loops over zonal faces.

      do J=js-1,je ; do i=is,ie ; if (G%vmask(i,J) > 0.5) then
        ! Determine the fluxes through the meridional faces.

        ! Find the acceptable range of tracer concentration around this face.
        if (nPv(i,J) >= 1) then
          Tr_min_face = min(Tr(m)%t(i,j,1), Tr(m)%t(i,j+1,1))
          Tr_max_face = max(Tr(m)%t(i,j,1), Tr(m)%t(i,j+1,1))
          do k=2,nkmb
            Tr_min_face = min(Tr_min_face, Tr(m)%t(i,j,k), Tr(m)%t(i,j+1,k))
            Tr_max_face = max(Tr_max_face, Tr(m)%t(i,j,k), Tr(m)%t(i,j+1,k))
          enddo

          ! Include the next two layers denser than the densest buffer layer.
          kLa = nkmb+1 ; if (max_kRho(i,j) < nz+1) kLa = max_kRho(i,j)
          kLb = kLa ; if (max_kRho(i,j) < nz) kLb = max_kRho(i,j)+1
          kRa = nkmb+1 ; if (max_kRho(i,j+1) < nz+1) kRa = max_kRho(i,j+1)
          kRb = kRa ; if (max_kRho(i,j+1) < nz) kRb = max_kRho(i,j+1)+1
          Tr_La = Tr_min_face ; Tr_Lb = Tr_La ; Tr_Ra = Tr_La ; Tr_Rb = Tr_La
          if (h(i,j,kLa) > h_exclude) Tr_La = Tr(m)%t(i,j,kLa)
          if (h(i,j,kLb) > h_exclude) Tr_La = Tr(m)%t(i,j,kLb)
          if (h(i,j+1,kRa) > h_exclude) Tr_Ra = Tr(m)%t(i,j+1,kRa)
          if (h(i,j+1,kRb) > h_exclude) Tr_Rb = Tr(m)%t(i,j+1,kRb)
          Tr_min_face = min(Tr_min_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
          Tr_max_face = max(Tr_max_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)

          ! Include all points in diffusive pairings at this face.
          do k=1,nPv(i,J)
            Tr_Lb = Tr(m)%t(i,j,k0b_Lv(J)%p(i,k)) ; Tr_Rb = Tr(m)%t(i,j+1,k0b_Rv(J)%p(i,k))
            Tr_La = Tr_Lb ; Tr_Ra = Tr_Rb
            if (deep_wt_Lv(J)%p(i,k) < 1.0) Tr_La = Tr(m)%t(i,j,k0a_Lv(J)%p(i,k))
            if (deep_wt_Rv(J)%p(i,k) < 1.0) Tr_Ra = Tr(m)%t(i,j+1,k0a_Rv(J)%p(i,k))
            Tr_min_face = min(Tr_min_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
            Tr_max_face = max(Tr_max_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
          enddo
        endif

        do k=1,nPv(i,J)
          kLb = k0b_Lv(J)%p(i,k) ; Tr_Lb = Tr(m)%t(i,j,kLb) ; Tr_av_L = Tr_Lb
          if (deep_wt_Lv(J)%p(i,k) < 1.0) then
            kLa = k0a_Lv(J)%p(i,k) ; Tr_La = Tr(m)%t(i,j,kLa)
            wt_b = deep_wt_Lv(J)%p(i,k)
            Tr_av_L = wt_b * Tr_Lb + (1.0-wt_b) * Tr_La
          endif

          kRb = k0b_Rv(J)%p(i,k) ; Tr_Rb = Tr(m)%t(i,j+1,kRb) ; Tr_av_R = Tr_Rb
          if (deep_wt_Rv(J)%p(i,k) < 1.0) then
            kRa = k0a_Rv(J)%p(i,k) ; Tr_Ra = Tr(m)%t(i,j+1,kRa)
            wt_b = deep_wt_Rv(J)%p(i,k)
            Tr_av_R = wt_b * Tr_Rb + (1.0-wt_b) * Tr_Ra
          endif

          h_L = hP_Lv(J)%p(i,k) ; h_R = hP_Rv(J)%p(i,k)
          Tr_flux = I_maxitt * ((2.0 * h_L * h_R) / (h_L + h_R)) * &
                    khdt_epi_y(i,J) * (Tr_av_L - Tr_av_R)
                    
          if (deep_wt_Lv(J)%p(i,k) >= 1.0) then
            tr_flux_conv(i,j,kLb) = tr_flux_conv(i,j,kLb) - Tr_flux
          else
            Tr_adj_vert = 0.0
            wt_b = deep_wt_Lv(J)%p(i,k) ; wt_a = 1.0 - wt_b
            vol = hP_Lv(J)%p(i,k) * G%DXDYh(i,j)

            !   Ensure that the tracer flux does not drive the tracer values
            ! outside of the range Tr_min_face <= Tr <= Tr_max_face.
            if (Tr_flux > 0.0) then
              if (Tr_La < Tr_Lb) then ; if (vol * (Tr_La-Tr_min_face) < Tr_flux) &
                Tr_adj_vert = -wt_a * min(Tr_flux - vol * (Tr_La-Tr_min_face), &
                                          (vol*wt_b) * (Tr_Lb - Tr_La))
              else ; if (vol*(Tr_Lb-Tr_min_face) < Tr_flux) &
                Tr_adj_vert = wt_b * min(Tr_flux - vol * (Tr_Lb-Tr_min_face), &
                                         (vol*wt_a) * (Tr_La - Tr_Lb))
              endif
            elseif (Tr_flux < 0.0) then
              if (Tr_La > Tr_Lb) then ; if (vol * (Tr_max_face-Tr_La) < -Tr_flux) &
                Tr_adj_vert = wt_a * min(-Tr_flux - vol * (Tr_max_face-Tr_La), &
                                         (vol*wt_b) * (Tr_La - Tr_Lb))
              else ; if (vol*(Tr_max_face-Tr_Lb) < -Tr_flux) &
                Tr_adj_vert = -wt_b * min(-Tr_flux - vol * (Tr_max_face-Tr_Lb), &
                                          (vol*wt_a)*(Tr_Lb - Tr_La))
              endif
            endif

            tr_flux_conv(i,j,kLa) = tr_flux_conv(i,j,kLa) - (wt_a*Tr_flux + Tr_adj_vert)
            tr_flux_conv(i,j,kLb) = tr_flux_conv(i,j,kLb) - (wt_b*Tr_flux - Tr_adj_vert)
          endif

          if (deep_wt_Rv(J)%p(i,k) >= 1.0) then
            tr_flux_conv(i,j+1,kRb) = tr_flux_conv(i,j+1,kRb) + Tr_flux
          else
            Tr_adj_vert = 0.0
            wt_b = deep_wt_Rv(J)%p(i,k) ; wt_a = 1.0 - wt_b
            vol = hP_Rv(J)%p(i,k) * G%DXDYh(i,j+1)

            !   Ensure that the tracer flux does not drive the tracer values
            ! outside of the range Tr_min_face <= Tr <= Tr_max_face.
            if (Tr_flux < 0.0) then
              if (Tr_Ra < Tr_Rb) then ; if (vol * (Tr_Ra-Tr_min_face) < -Tr_flux) &
                Tr_adj_vert = -wt_a * min(-Tr_flux - vol * (Tr_Ra-Tr_min_face), &
                                          (vol*wt_b) * (Tr_Rb - Tr_Ra))
              else ; if (vol*(Tr_Rb-Tr_min_face) < (-Tr_flux)) &
                Tr_adj_vert = wt_b * min(-Tr_flux - vol * (Tr_Rb-Tr_min_face), &
                                         (vol*wt_a) * (Tr_Ra - Tr_Rb))
              endif
            elseif (Tr_flux > 0.0) then
              if (Tr_Ra > Tr_Rb) then ; if (vol * (Tr_max_face-Tr_Ra) < Tr_flux) &
                Tr_adj_vert = wt_a * min(Tr_flux - vol * (Tr_max_face-Tr_Ra), &
                                         (vol*wt_b) * (Tr_Ra - Tr_Rb))
              else ; if (vol*(Tr_max_face-Tr_Rb) < Tr_flux) &
                Tr_adj_vert = -wt_b * min(Tr_flux - vol * (Tr_max_face-Tr_Rb), &
                                          (vol*wt_a)*(Tr_Rb - Tr_Ra))
              endif
            endif

            tr_flux_conv(i,j+1,kRa) = tr_flux_conv(i,j+1,kRa) + &
                                            (wt_a*Tr_flux - Tr_adj_vert)
            tr_flux_conv(i,j+1,kRb) = tr_flux_conv(i,j+1,kRb) + &
                                            (wt_b*Tr_flux + Tr_adj_vert)
          endif
          if (associated(Tr(m)%df2d_y)) &
            Tr(m)%df2d_y(i,J) = Tr(m)%df2d_y(i,J) + Tr_flux * Idt
        enddo ! Loop over pairings at faces.
      endif ; enddo ; enddo ! i- & j- loops over meridional faces.


      do k=1,PEmax_kRho ; do j=js,je ; do i=is,ie
        if ((G%hmask(i,j) > 0.5) .and. (h(i,j,k) > 0.0)) then
          Tr(m)%t(i,j,k) = Tr(m)%t(i,j,k) + tr_flux_conv(i,j,k) / &
                                            (h(i,j,k)*G%DXDYh(i,j))
          tr_flux_conv(i,j,k) = 0.0
        endif
      enddo ; enddo ; enddo

    enddo ! Loop over tracers
  enddo ! Loop over iterations

  do j=js,je
    deallocate(deep_wt_Lu(j)%p) ; deallocate(deep_wt_Ru(j)%p)
    deallocate(Hp_Lu(j)%p)  ; deallocate(Hp_Ru(j)%p)
    deallocate(k0a_Lu(j)%p) ; deallocate(k0a_Ru(j)%p)
    deallocate(k0b_Lu(j)%p) ; deallocate(k0b_Ru(j)%p)
  enddo

  do J=js-1,je
    deallocate(deep_wt_Lv(J)%p) ; deallocate(deep_wt_Rv(J)%p)
    deallocate(Hp_Lv(J)%p)  ; deallocate(Hp_Rv(J)%p)
    deallocate(k0a_Lv(J)%p) ; deallocate(k0a_Rv(J)%p)
    deallocate(k0b_Lv(J)%p) ; deallocate(k0b_Rv(J)%p)
  enddo

end subroutine tracer_epipycnal_ML_diff

subroutine tracer_vertdiff(h_old, ea, eb, dt, tr, G, &
                           sfc_flux, btm_flux, btm_reservoir, sink_rate)
  real, dimension(NXMEM_,NYMEM_,NKMEM_), intent(in)    :: h_old, ea, eb
  real, dimension(NXMEM_,NYMEM_,NKMEM_), intent(inout) :: tr
  real,                                  intent(in)    :: dt
  type(ocean_grid_type),                 intent(in)    :: G
  real, dimension(NXMEM_,NYMEM_), optional, intent(in) :: sfc_flux
  real, dimension(NXMEM_,NYMEM_), optional, intent(in) :: btm_flux
  real, dimension(NXMEM_,NYMEM_), optional, intent(inout) :: btm_reservoir
  real,                           optional, intent(in) :: sink_rate
! Arguments: h_old -  Layer thickness before entrainment, in m or kg m-2.
!  (in)      ea - The amount of fluid entrained from the layer above, in the
!                 same units as h_old, i.e. m or kg m-2.
!  (in)      eb - The amount of fluid entrained from the layer below, in the
!                 same units as h_old, i.e. m or kg m-2
!  (inout)   tr - The tracer concentration, in concentration units (CU).
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!  (in,opt)  sfc_flux - The surface flux of the tracer, in CU kg m-2 s-1.
!  (in,opt)  btm_flux - The (negative upward) bottom flux of the tracer,
!                       in units of CU kg m-2 s-1.
!  (inout,opt) btm_reservoir - The amount of tracer in a bottom reservoir, in
!                              units of CU kg m-2. (was CU m)
!  (in,opt)  sink_rate - The rate at which the tracer sinks, in m s-1.

!   This subroutine solves a tridiagonal equation for the final tracer
! concentrations after the dual-entrainments, and possibly sinking or surface
! and bottom sources, are applied.  The sinking is implemented with an
! fully implicit upwind advection scheme.
 
  real :: sink_dist ! The distance the tracer sinks in a time step, in m or kg m-2.
  real, dimension(SZI_(G)) :: &
    sfc_src, &      ! The time-integrated surface source of the tracer, in
                    ! units of m or kg m-2 times a concentration.
    btm_src, &      ! The time-integrated bottom source of the tracer, in
                    ! units of m or kg m-2  times a concentration.
    b1, &           ! b1 is used by the tridiagonal solver, in m-1 or m2 kg-1.
    d1              ! d1=1-c1 is used by the tridiagonal solver, nondimensional.
  real :: c1(SZI_(G),SZK_(G))     ! c1 is used by the tridiagonal solver, ND.
  real :: h_minus_dsink(SZI_(G),SZK_(G))  ! The layer thickness minus the
                    ! difference in sinking rates across the layer, in m or kg m-2.
                    ! By construction, 0 <= h_minus_dsink < h_old.
  real :: sink(SZI_(G),SZK_(G)+1) ! The tracer's sinking distances at the
                    ! interfaces, limited to prevent characteristics from
                    ! crossing within a single timestep, in m or kg m-2.
  real :: b_denom_1 ! The first term in the denominator of b1, in m or kg m-2.
  real :: h_tr      ! h_tr is h at tracer points with a h_neglect added to
                    ! ensure positive definiteness, in m or kg m-2.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in m.
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  h_neglect = G%H_subroundoff

  sink_dist = 0.0
  if (present(sink_rate)) sink_dist = (dt*sink_rate) * G%m_to_H
  do i=is,ie ; sfc_src(i) = 0.0 ; btm_src(i) = 0.0 ; enddo

  if (present(sink_rate)) then
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
      do i=is,ie ; if (G%hmask(i,j) > 0.5) then
        b_denom_1 = h_minus_dsink(i,1) + ea(i,j,1) + h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,1))
        d1(i) = b_denom_1 * b1(i)
        if (present(sfc_flux)) sfc_src(i) = (sfc_flux(i,j)*dt) * G%kg_m2_to_H
        h_tr = h_old(i,j,1) + h_neglect
        tr(i,j,1) = b1(i)*(h_tr*tr(i,j,1) + sfc_src(i))
      endif ; enddo
      do k=2,nz-1 ; do i=is,ie ; if (G%hmask(i,j) > 0.5) then
        c1(i,k) = eb(i,j,k-1) * b1(i)
        b_denom_1 = h_minus_dsink(i,k) + d1(i) * (ea(i,j,k) + sink(i,K)) + &
                    h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
        d1(i) = b_denom_1 * b1(i)
        h_tr = h_old(i,j,k) + h_neglect
        tr(i,j,k) = b1(i) * (h_tr * tr(i,j,k) + &
                             (ea(i,j,k) + sink(i,K)) * tr(i,j,k-1))
      endif ; enddo ; enddo
      do i=is,ie ; if (G%hmask(i,j) > 0.5) then
        c1(i,nz) = eb(i,j,nz-1) * b1(i)
        b_denom_1 = h_minus_dsink(i,nz) + d1(i) * (ea(i,j,nz) + sink(i,nz)) + &
                    h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,nz))
        h_tr = h_old(i,j,nz) + h_neglect
        if (present(btm_flux)) btm_src(i) = (-btm_flux(i,j)*dt) * G%kg_m2_to_H
        tr(i,j,nz) = b1(i) * ((h_tr * tr(i,j,nz) + btm_src(i)) + &
                              (ea(i,j,nz) + sink(i,nz)) * tr(i,j,nz-1))
      endif ; enddo
      if (present(btm_reservoir)) then ; do i=is,ie ; if (G%hmask(i,j)>0.5) then
        btm_reservoir(i,j) = btm_reservoir(i,j) + &
                             (sink(i,nz+1)*tr(i,j,nz)) * G%H_to_kg_m2
      endif ; enddo ; endif

      do k=nz-1,1,-1 ; do i=is,ie ; if (G%hmask(i,j) > 0.5) then
        tr(i,j,k) = tr(i,j,k) + c1(i,k+1)*tr(i,j,k+1)
      endif ; enddo ; enddo
    enddo
  else
    do j=js,je
      do i=is,ie ; if (G%hmask(i,j) > 0.5) then
        h_tr = h_old(i,j,1) + h_neglect
        b_denom_1 = h_tr + ea(i,j,1)
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,1))
        d1(i) = b_denom_1 * b1(i)
        if (present(sfc_flux)) sfc_src(i) = (sfc_flux(i,j)*dt) * G%kg_m2_to_H
        tr(i,j,1) = b1(i)*(h_tr*tr(i,j,1) + sfc_src(i))
      endif ; enddo
      do k=2,nz-1 ; do i=is,ie ; if (G%hmask(i,j) > 0.5) then
        c1(i,k) = eb(i,j,k-1) * b1(i)
        h_tr = h_old(i,j,k) + h_neglect
        b_denom_1 = h_tr + d1(i) * ea(i,j,k)
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
        d1(i) = b_denom_1 * b1(i)
        tr(i,j,k) = b1(i) * (h_tr * tr(i,j,k) + ea(i,j,k) * tr(i,j,k-1))
      endif ; enddo ; enddo
      do i=is,ie ; if (G%hmask(i,j) > 0.5) then
        c1(i,nz) = eb(i,j,nz-1) * b1(i)
        h_tr = h_old(i,j,nz) + h_neglect
        b1(i) = 1.0 / (h_tr + d1(i) * ea(i,j,nz) + eb(i,j,nz))
        if (present(btm_flux)) btm_src(i) = (-btm_flux(i,j)*dt) * G%kg_m2_to_H
        tr(i,j,nz) = b1(i) * ((h_tr * tr(i,j,nz) + btm_src(i)) + &
                              ea(i,j,nz) * tr(i,j,nz-1))
      endif ; enddo
      do k=nz-1,1,-1 ; do i=is,ie ; if (G%hmask(i,j) > 0.5) then
        tr(i,j,k) = tr(i,j,k) + c1(i,k+1)*tr(i,j,k+1)
      endif ; enddo ; enddo
    enddo
  endif

end subroutine tracer_vertdiff

subroutine MOM_tracer_chksum(mesg, Tr, ntr, G)
  character(len=*),         intent(in) :: mesg
  type(tracer),             intent(in) :: Tr(:)
  integer,                  intent(in) :: ntr
  type(ocean_grid_type),    intent(in) :: G
!   This subroutine writes out chksums for the model's thermodynamic state
! variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      Tr - An array of all of the registered tracers.
!  (in)      ntr - The number of registered tracers.
!  (in)      G - The ocean's grid structure.
  integer :: is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  do m=1,ntr
    call hchksum(Tr(m)%t, mesg//trim(Tr(m)%name), G)
  enddo
end subroutine MOM_tracer_chksum

subroutine advect_tracer_init(param_file, CS)
  type(param_file_type), intent(in)  :: param_file
  type(advect_tracer_CS), pointer    :: CS
! Arguments: param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
  integer, save :: init_calls = 0
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'
  character(len=40)  :: mod = "MOM_tracer" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (.not.associated(CS)) then ; allocate(CS)
  else ; return ; endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, tagname, "")
  call get_param(param_file, mod, "KHTR", CS%KhTr, &
                 "The background along-isopycnal tracer diffusivity.", &
                 units="m2 s-1", default=0.0)
  call get_param(param_file, mod, "KHTR_SLOPE_CFF", CS%KhTr_Slope_Cff, &
                 "The scaling coefficient for along-isopycnal tracer \n"//&
                 "diffusivity using a shear-based (Visbeck-like) \n"//&
                 "parameterization.  A non-zero value enables this param.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mod, "KHTR_MIN", CS%KhTr_Min, &
                 "The minimum along-isopycnal tracer diffusivity.", &
                 units="m2 s-1", default=0.0)
  call get_param(param_file, mod, "KHTR_MAX", CS%KhTr_Max, &
                 "The maximum along-isopycnal tracer diffusivity.", &
                 units="m2 s-1", default=0.0)
  call get_param(param_file, mod, "KHTR_PASSIVITY_COEFF", CS%KhTr_passivity_coeff, &
               "The coefficient that scales deformation radius over \n"//&
               "grid-spacing in passivity, where passiviity is the ratio \n"//&
               "between along isopycnal mxiing of tracers to thickness mixing. \n"//&
               "A non-zero value enables this parameterization.", &
               units="nondim", default=0.0)
  call get_param(param_file, mod, "KHTR_PASSIVITY_MIN", CS%KhTr_passivity_min, &
               "The minimum passivity which is the ratio between \n"//&
               "along isopycnal mxiing of tracers to thickness mixing. \n", &
               units="nondim", default=0.5)
  call get_param(param_file, mod, "DT", CS%dt, fail_if_missing=.true., &
          desc="The (baroclinic) dynamics time step.", units="s")


  call get_param(param_file, mod, "DIFFUSE_ML_TO_INTERIOR", CS%Diffuse_ML_interior, &
                 "If true, enable epipycnal mixing between the surface \n"//&
                 "boundary layer and the interior.", default=.false.)
  call get_param(param_file, mod, "CHECK_DIFFUSIVE_CFL", CS%check_diffusive_CFL, &
                 "If true, use enough iterations the diffusion to ensure \n"//&
                 "that the diffusive equivalent of the CFL limit is not \n"//&
                 "violated.  If false, always use 1 iteration.", default=.false.)
  CS%ML_KhTR_scale = 1.0
  if (CS%Diffuse_ML_interior) then
    call get_param(param_file, mod, "ML_KHTR_SCALE", CS%ML_KhTR_scale, &
                 "With Diffuse_ML_interior, the ratio of the truly \n"//&
                 "horizontal diffusivity in the mixed layer to the \n"//&
                 "epipycnal diffusivity.  The valid range is 0 to 1.", &
                 units="nondim", default=1.0)
  endif

  call get_param(param_file, mod, "DEBUG", CS%debug, default=.false.)

  init_calls = init_calls + 1
  if (init_calls == 1) then

    id_clock_advect = cpu_clock_id('(Ocean advect tracer)', grain=CLOCK_MODULE)
    id_clock_diffuse = cpu_clock_id('(Ocean diffuse tracer)', grain=CLOCK_MODULE)
    id_clock_epimix = cpu_clock_id('(Ocean epipycnal diffuse tracer)', grain=CLOCK_MODULE)
    id_clock_pass = cpu_clock_id('(Ocean tracer halo updates)', grain=CLOCK_ROUTINE)
    id_clock_sync = cpu_clock_id('(Ocean tracer global synch)', grain=CLOCK_ROUTINE)
  else
    write(mesg,'("Advect_tracer_init called ",I3, &
      &" times with different control structure pointers.")') init_calls
    if (is_root_pe()) call MOM_error(WARNING,"MOM_tracer"//mesg)
  endif

  CS%id_KhTr_u = -1 ; CS%id_KhTr_v = -1

end subroutine advect_tracer_init

subroutine advect_tracer_diag_init(Time, G, diag, CS)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(diag_ptrs), target, intent(inout) :: diag
  type(advect_tracer_CS),  pointer       :: CS
!   This subroutine sets up any diagnostics that are common to all tracers.  It
! needs to be called separately from advect_tracer_init because the time was not
! known yet hen advect_tracer_init was called.
!
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  CS - A pointer to the control structure for this module

  !   If the control structure is not associated, there are probably no tracers,
  ! so return without giving an error.
  if (.not.associated(CS)) return

  CS%diag => diag
  CS%id_KhTr_u = register_diag_field('ocean_model', 'KhTr_u', G%axesu1, Time, &
     'Epipycnal tracer diffusivity at zonal faces of tracer cell', 'meter2 second-1')
  CS%id_KhTr_v = register_diag_field('ocean_model', 'KhTr_v', G%axesv1, Time, &
     'Epipycnal tracer diffusivity at meridional faces of tracer cell', 'meter2 second-1')

end subroutine advect_tracer_diag_init

subroutine tracer_end(CS)
  type(advect_tracer_CS), pointer :: CS
  if (associated(CS)) deallocate(CS)
end subroutine tracer_end

end module MOM_tracer
