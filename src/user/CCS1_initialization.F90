module CCS1_initialization
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

!***********************************************************************
!*                                                                     *
!*  The module configures the model for the "CCS1" experiment.         *
!*  CCS1 = California Current System (one of many)                     *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type, add_tracer_OBC_values
use MOM_variables, only : thermo_var_ptrs, ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

implicit none ; private

#include <MOM_memory.h>

public CCS1_initialize_sponges
public CCS1_set_Open_Bdry_Conds
public CCS1_set_Open_Bdry_Vals

contains

! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine CCS1_initialize_sponges(G, GV, tv, PF, CSp)
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  type(thermo_var_ptrs), intent(in) :: tv
  type(param_file_type), intent(in) :: PF
  type(sponge_CS),       pointer    :: CSp
!   This subroutine sets the inverse restoration time (Idamp), and   !
! the values towards which the interface heights and an arbitrary    !
! number of tracers should be restored within each sponge. The       !
! interface height is always subject to damping, and must always be  !
! the first registered field.                                        !

! Arguments: G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in)      PF - A structure indicating the open file to parse for
!                 model parameter values.
!  (in/out)  CSp - A pointer that is set to point to the control structure
!                  for this module

  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! A temporary array for eta.
  real :: temp(SZI_(G),SZJ_(G),SZK_(G))  ! A temporary array for other variables. !
  real :: Idamp(SZI_(G),SZJ_(G))    ! The inverse damping rate, in s-1.

  real :: H0(SZK_(G))
  real :: min_depth
  real :: damp, e_dense, damp_new
  character(len=40)  :: mod = "CCS1_initialize_sponges" ! This subroutine's name.
  integer :: i, j, k, isc, iec, jsc, jec, isd, ied, jsd, jed, nz

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  eta(:,:,:) = 0.0 ; temp(:,:,:) = 0.0 ; Idamp(:,:) = 0.0

!  Here the inverse damping time, in s-1, is set. Set Idamp to 0     !
!  wherever there is no sponge, and the subroutines that are called  !
!  will automatically set up the sponges only where Idamp is positive!
!  and mask2dT is 1.                                                   !

!   Set up sponges for CCS1 configuration
  call get_param(PF, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)

  H0(1) = 0.0
  do k=2,nz ; H0(k) = -(real(k-1)-0.5)*G%max_depth/real(nz-1) ; enddo
  do i=isc,iec; do j=jsc,jec
    if (G%geoLonT(i,j) < 100.0) then ; damp = 10.0
    elseif (G%geoLonT(i,j) < 200.0) then
      damp = 10.0*(200.0-G%geoLonT(i,j))/100.0
    else ; damp=0.0
    endif

    if (G%geoLonT(i,j) > 1400.0) then ; damp_new = 10.0
    elseif (G%geoLonT(i,j) > 1300.0) then
       damp_new = 10.0*(G%geoLonT(i,j)-1300.0)/100.0
    else ; damp_new = 0.0
    endif

    if (damp <= damp_new) damp=damp_new

    ! These will be stretched inside of apply_sponge, so they can be in
    ! depth space for Boussinesq or non-Boussinesq models.
    eta(i,j,1) = 0.0
    do k=2,nz
!     eta(i,j,K)=max(H0(k), -G%bathyT(i,j), GV%Angstrom_z*(nz-k+1)-G%bathyT(i,j))
      e_dense = -G%bathyT(i,j)
      if (e_dense >= H0(k)) then ; eta(i,j,K) = e_dense
      else ; eta(i,j,K) = H0(k) ; endif
      if (eta(i,j,K) < GV%Angstrom_z*(nz-k+1)-G%bathyT(i,j)) &
          eta(i,j,K) = GV%Angstrom_z*(nz-k+1)-G%bathyT(i,j)
    enddo
    eta(i,j,nz+1) = -G%bathyT(i,j)

    if (G%bathyT(i,j) > min_depth) then
      Idamp(i,j) = damp/86400.0
    else ; Idamp(i,j) = 0.0 ; endif
  enddo ; enddo

!  This call sets up the damping rates and interface heights.
!  This sets the inverse damping timescale fields in the sponges.    !
  call initialize_sponge(Idamp, eta, G, PF, CSp)

!   Now register all of the fields which are damped in the sponge.   !
! By default, momentum is advected vertically within the sponge, but !
! momentum is typically not damped within the sponge.                !

! At this point, the CCS1 configuration is done. The following are here as a
! template for other configurations.

!  The remaining calls to set_up_sponge_field can be in any order. !
  if ( associated(tv%T) ) then
    call MOM_error(FATAL,"CCS1_initialize_sponges is not set up for use with"//&
                         " a temperature defined.")
    ! This should use the target values of T in temp.
    call set_up_sponge_field(temp, tv%T, G, nz, CSp)
    ! This should use the target values of S in temp.
    call set_up_sponge_field(temp, tv%S, G, nz, CSp)
  endif

end subroutine CCS1_initialize_sponges
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine CCS1_set_Open_Bdry_Conds(OBC, tv, G, GV, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC
  type(thermo_var_ptrs),      intent(in) :: tv
  type(ocean_grid_type),      intent(in) :: G
  type(verticalGrid_type),    intent(in) :: GV
  type(param_file_type),      intent(in) :: param_file
  type(tracer_registry_type), pointer    :: tr_Reg
!   This subroutine sets the structure for open boundary conditions.
! This particular example is for the CCS1, or any rectangular domain with open
! boundaries at the edges.

! Arguments: OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

  logical :: any_OBC        ! Set to true if any points in this subdomain use
                            ! open boundary conditions.
  logical, pointer, dimension(:,:) :: &
    OBC_mask_u => NULL(), & ! These arrays are true at zonal or meridional
    OBC_mask_v => NULL()    ! velocity points that have prescribed open boundary
                            ! conditions.
  real, pointer, dimension(:,:,:) :: &
    OBC_T_u => NULL(), &    ! These arrays should be allocated and set to
    OBC_T_v => NULL(), &    ! specify the values of T and S that should come
    OBC_S_u => NULL(), &    ! in through u- and v- points through the open
    OBC_S_v => NULL()       ! boundary conditions, in C and psu.
  logical :: apply_OBC_u, apply_OBC_v
  ! The following variables are used to set the target temperature and salinity.
  real :: T0(SZK_(G)), S0(SZK_(G))
  real :: pres(SZK_(G))      ! An array of the reference pressure in Pa.
  real :: drho_dT(SZK_(G))   ! Derivative of density with temperature in kg m-3 K-1.                              !
  real :: drho_dS(SZK_(G))   ! Derivative of density with salinity in kg m-3 PSU-1.                             !
  real :: rho_guess(SZK_(G)) ! Potential density at T0 & S0 in kg m-3.
  ! The following variables are used to set up the transport in the CCS1 example.
  real :: tr_0, y1, y2, tr_k, rst, rsb, rc, v_k, lon_im1
  character(len=40)  :: mod = "CCS1_set_Open_Bdry_Conds" ! This subroutine's name.
  integer :: i, j, k, itt, isc, iec, jsc, jec, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB, IscB, IecB, JscB, JecB

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  IscB = G%IscB ; IecB = G%IecB ; JscB = G%JscB ; JecB = G%JecB

  call get_param(param_file, mod, "APPLY_OBC_U", apply_OBC_u, &
                 "If true, open boundary conditions may be set at some \n"//&
                 "u-points, with the configuration controlled by OBC_CONFIG", &
                 default=.false.)
  call get_param(param_file, mod, "APPLY_OBC_V", apply_OBC_v, &
                 "If true, open boundary conditions may be set at some \n"//&
                 "v-points, with the configuration controlled by OBC_CONFIG", &
                 default=.false.)

  if (apply_OBC_u) then
    ! Determine where u points are applied.
    allocate(OBC_mask_u(IsdB:IedB,jsd:jed)) ; OBC_mask_u(:,:) = .false.
    any_OBC = .false.
    ! Check for edges of full domain
    ! West
!    if (is_root_pe()) print *, 'inside CCS1_set_open_bdry_conds 1', &
!      apply_OBC_u, apply_OBC_v, G%isd_global, isd, isc, iec, IscB, IecB, G%mask2dCu(iscB,jsc)
    if (G%isd_global - isd == isd - isc) then
      do j=jsc,jec
        if (G%mask2dT(isc,j) == 1.0) then
          OBC_mask_u(isc-1,j) = .true. ; any_OBC = .true.
        endif
      enddo
    endif
    ! East - not wanted for CCS1
!    if (G%isd_global - isd + iec == G%Domain%niglobal) then
!      do j=jsc,jec
!        if (G%mask2dT(Iec,j) == 1.0) then
!          OBC_mask_u(Iec,j) = .true. ; any_OBC = .true.
!        endif
!      enddo
!    endif
    ! South
    if (G%jsd_global - jsd == jsd - jsc) then
      do I=IscB,IecB
        if (G%mask2dCu(I,jsc) == 1.0) then
          OBC_mask_u(I,jsc-1) = .true. ; any_OBC = .true.
        endif
      enddo
    endif
    ! North
    if (G%jsd_global - jsd + jec == G%Domain%njglobal) then
      do I=IscB,IecB
        if (G%mask2dCu(I,jec) == 1.0) then
          OBC_mask_u(I,jec+1) = .true. ; any_OBC = .true.
        endif
      enddo
    endif
    print *, 'inside CCS1_set_open_bdry_conds ', any_OBC, G%jsd_global, jsd, &
        jsc, jec, jscB, jecB, G%Domain%njglobal
    if (.not.any_OBC) then
      ! This processor has no u points at which open boundary conditions are
      ! to be applied.
      apply_OBC_u = .false.
      deallocate(OBC_mask_u)
    endif
  endif
  if (apply_OBC_v) then
    ! Determine where v points are applied.
    allocate(OBC_mask_v(isd:ied,JsdB:JedB)) ; OBC_mask_v(:,:) = .false.
    any_OBC = .false.
    ! Check for edges of full domain
    ! West
    if (G%isd_global - isd == isd - isc) then
      do J=JscB,JecB
        if (G%mask2dCv(isc,J) == 1) then
          OBC_mask_v(isc-1,J) = .true. ; any_OBC = .true.
        endif
      enddo
    endif
    ! East - not for CCS
!    if (G%isd_global - isd + iec == G%Domain%niglobal) then
!      do J=JscB,JecB
!        if (G%mask2dCv(iec,J) == 1) then
!          OBC_mask_v(iec+1,J) = .true. ; any_OBC = .true.
!        endif
!      enddo
!    endif
    ! South
    if (G%jsd_global - jsd == jsd - jsc) then
      do i=isc,iec
        if (G%mask2dT(i,jsc) == 1) then
          OBC_mask_v(i,jsc-1) = .true. ; any_OBC = .true.
        endif
      enddo
    endif
    ! North
    if (G%jsd_global - jsd + jec == G%Domain%njglobal) then
      do i=isc,iec
        if (G%mask2dT(i,jec) == 1) then
          OBC_mask_v(i,Jec) = .true. ; any_OBC = .true.
        endif
      enddo
    endif
    if (.not.any_OBC) then
      ! This processor has no v points at which open boundary conditions are
      ! to be applied.
      apply_OBC_v = .false.
      deallocate(OBC_mask_v)
    endif
  endif

  if (is_root_pe()) print *, 'inside CCS1_set_open_bdry_conds', apply_OBC_u, apply_OBC_v
  if (.not.(apply_OBC_u .or. apply_OBC_v)) return

  if (.not.associated(OBC)) allocate(OBC)   

  if (apply_OBC_u) then
    OBC%apply_OBC_u = .true.
    OBC%OBC_mask_u => OBC_mask_u
    allocate(OBC%u(IsdB:IedB,jsd:jed,nz)) ; OBC%u(:,:,:) = 0.0
    allocate(OBC%uh(IsdB:IedB,jsd:jed,nz)) ; OBC%uh(:,:,:) = 0.0
    allocate(OBC%OBC_kind_u(IsdB:IedB,jsd:jed)) ; OBC%OBC_kind_u(:,:) = OBC_NONE
    do j=jsd,jed ; do I=IsdB,IedB
      if (OBC%OBC_mask_u(I,j)) OBC%OBC_kind_u(I,j) = OBC_SIMPLE
    enddo ; enddo
  endif
  if (apply_OBC_v) then
    OBC%apply_OBC_v = .true.
    OBC%OBC_mask_v => OBC_mask_v
    allocate(OBC%v(isd:ied,JsdB:JedB,nz)) ; OBC%v(:,:,:) = 0.0
    allocate(OBC%vh(isd:ied,JsdB:JedB,nz)) ; OBC%vh(:,:,:) = 0.0
    allocate(OBC%OBC_kind_v(isd:ied,JsdB:JedB)) ; OBC%OBC_kind_v(:,:) = OBC_NONE
    do J=JsdB,JedB ; do i=isd,ied
      if (OBC%OBC_mask_v(i,J)) OBC%OBC_kind_v(i,J) = OBC_SIMPLE
    enddo ; enddo
  endif

  !   The inflow values of temperature and salinity also need to be set here if
  ! these variables are used.  The following code is just a naive example.
  if (apply_OBC_u .or. apply_OBC_v) then
    if (associated(tv%S)) then
      ! In this example, all S inflows have values of 35 psu.
      call add_tracer_OBC_values("S", tr_Reg, OBC_inflow=35.0)
    endif
    if (associated(tv%T)) then
      ! In this example, the T values are set to be consistent with the layer
      ! target density and a salinity of 35 psu.  This code is taken from
      ! USER_initialize_temp_sal.
      pres(:) = tv%P_Ref ; S0(:) = 35.0 ; T0(1) = 25.0
      call calculate_density(T0(1),S0(1),pres(1),rho_guess(1),tv%eqn_of_state)
      call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,1,tv%eqn_of_state)

      do k=1,nz ; T0(k) = T0(1) + (GV%Rlay(k)-rho_guess(1)) / drho_dT(1) ; enddo
      do itt=1,6
        call calculate_density(T0,S0,pres,rho_guess,1,nz,tv%eqn_of_state)
        call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,nz,tv%eqn_of_state)
        do k=1,nz ; T0(k) = T0(k) + (GV%Rlay(k)-rho_guess(k)) / drho_dT(k) ; enddo
      enddo

      if (apply_OBC_u) then
        allocate(OBC_T_u(IsdB:IedB,jsd:jed,nz))
        do k=1,nz ; do j=jsd,jed ; do I=IsdB,IedB
          OBC_T_u(I,j,k) = T0(k)
        enddo ; enddo ; enddo
      endif
      if (apply_OBC_v) then
        allocate(OBC_T_v(isd:ied,JsdB:JedB,nz))
        do k=1,nz ; do J=JsdB,JedB ; do i=isd,ied
          OBC_T_v(i,J,k) = T0(k)
        enddo ; enddo ; enddo
      endif
      call add_tracer_OBC_values("T", tr_Reg, OBC_in_u=OBC_T_u, &
                                              OBC_in_v=OBC_T_v)
    endif
  endif

end subroutine CCS1_set_Open_Bdry_Conds
! -----------------------------------------------------------------------------
subroutine CCS1_set_Open_Bdry_Vals(OBC, u, v, h, tv, G, GV, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in)   :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in)   :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)   :: h
  type(thermo_var_ptrs),      intent(in) :: tv
  type(ocean_grid_type),      intent(in) :: G
  type(verticalGrid_type),    intent(in) :: GV
  type(param_file_type),      intent(in) :: param_file
  type(tracer_registry_type), pointer    :: tr_Reg
!   This subroutine sets the properties of flow at open boundary conditions.
! This particular example applies the CCS1 initial conditions to the boundary
! for all time.

! Arguments: OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in)      u  - Zonal velocity, in m s-1.
!  (in)      v  - Meridional velocity, in m s-1.
!  (in)      h  - Layer thickness, in m.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

  logical :: any_OBC        ! Set to true if any points in this subdomain use
                            ! open boundary conditions.
  logical, pointer, dimension(:,:) :: &
    OBC_mask_u => NULL(), & ! These arrays are true at zonal or meridional
    OBC_mask_v => NULL()    ! velocity points that have prescribed open boundary
                            ! conditions.
  real, pointer, dimension(:,:,:) :: &
    OBC_T_u => NULL(), &    ! These arrays should be allocated and set to
    OBC_T_v => NULL(), &    ! specify the values of T and S that should come
    OBC_S_u => NULL(), &    ! in through u- and v- points through the open
    OBC_S_v => NULL()       ! boundary conditions, in C and psu.
  logical :: apply_OBC_u, apply_OBC_v
  ! The following variables are used to set the target temperature and salinity.
  real :: T0(SZK_(G)), S0(SZK_(G))
  real :: pres(SZK_(G))      ! An array of the reference pressure in Pa.
  real :: drho_dT(SZK_(G))   ! Derivative of density with temperature in kg m-3 K-1.                              !
  real :: drho_dS(SZK_(G))   ! Derivative of density with salinity in kg m-3 PSU-1.                             !
  real :: rho_guess(SZK_(G)) ! Potential density at T0 & S0 in kg m-3.
  ! The following variables are used to set up the transport in the CCS1 example.
  real :: tr_0, y1, y2, tr_k, rst, rsb, rc, v_k, lon_im1
  character(len=40)  :: mod = "CCS1_set_Open_Bdry_Conds" ! This subroutine's name.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) return
  if (.not.(OBC%apply_OBC_u .or. OBC%apply_OBC_v)) return

  if (OBC%apply_OBC_u) then
    do k=1,nz ; do j=jsd,jed ; do I=IsdB,IedB
      if (OBC%OBC_mask_u(I,j)) then
        ! An appropriate expression for the zonal inflow velocities and
        ! transports should go here.
        OBC%uh(I,j,k) = h(I,j,k) * u(I,j,k)
        OBC%u(I,j,k) = u(I,j,k)
      else
        OBC%uh(I,j,k) = 0.0 ; OBC%u(I,j,k) = 0.0
      endif
    enddo ; enddo ; enddo
  endif

  if (OBC%apply_OBC_v) then
    do k=1,nz
      do J=JsdB,JedB ; do i=isd,ied
        if (OBC%OBC_mask_v(i,J)) then
          OBC%vh(i,J,k) = v(i,J,k)*h(i,J,k)
          OBC%v(i,J,k) = v(i,J,k)
        else
          OBC%vh(i,J,k) = 0.0 ; OBC%v(i,J,k) = 0.0
        endif
      enddo ; enddo
    enddo
  endif

  !   The inflow values of temperature and salinity also need to be set here if
  ! these variables are used.  The following code is just a naive example.
  if (associated(tv%S)) then
    ! In this example, all S inflows have values of 35 psu.
!    call add_tracer_OBC_values("S", tr_Reg, OBC_inflow=35.0)
  endif
  if (associated(tv%T)) then
    ! In this example, the T values are set to be consistent with the layer
  endif

end subroutine CCS1_set_Open_Bdry_Vals
! -----------------------------------------------------------------------------

end module CCS1_initialization
