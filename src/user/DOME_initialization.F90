module DOME_initialization
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

use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_open_boundary, only : OBC_segment_type
use MOM_tracer_registry, only : tracer_registry_type, add_tracer_OBC_values
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

implicit none ; private

#include <MOM_memory.h>

public DOME_initialize_topography
public DOME_initialize_thickness
public DOME_initialize_sponges
public DOME_set_OBC_data

contains

! -----------------------------------------------------------------------------
!> This subroutine sets up the DOME topography
subroutine DOME_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m

  real :: min_depth ! The minimum and maximum depths in m.
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "DOME_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call MOM_mesg("  DOME_initialization.F90, DOME_initialize_topography: setting topography", 5)

  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)

  do j=js,je ; do i=is,ie
    if (G%geoLatT(i,j) < 600.0) then
      if (G%geoLatT(i,j) < 300.0) then
        D(i,j)=max_depth
      else
        D(i,j)=max_depth-10.0*(G%geoLatT(i,j)-300.0)
      endif
    else
      if ((G%geoLonT(i,j) > 1000.0).AND.(G%geoLonT(i,j) < 1100.0)) then
        D(i,j)=600.0
      else
        D(i,j)=0.5*min_depth
      endif
    endif

    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
  enddo ; enddo

end subroutine DOME_initialize_topography
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> This subroutine initializes layer thicknesses for the DOME experiment
subroutine DOME_initialize_thickness(h, G, GV, param_file)
  type(ocean_grid_type),   intent(in) :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV !< The ocean's vertical grid structure.
  real, intent(out), dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h !< The thickness that is being initialized.
  type(param_file_type),   intent(in) :: param_file !< A structure indicating the open file
                                                    !! to parse for model parameter values.

  real :: e0(SZK_(G)+1)     ! The resting interface heights, in m, usually !
                            ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)  ! Interface height relative to the sea surface !
                            ! positive upward, in m.                       !
  character(len=40)  :: mod = "DOME_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("  DOME_initialization.F90, DOME_initialize_thickness: setting thickness", 5)

  e0(1)=0.0
  do k=2,nz
    e0(K) = -G%max_depth * (real(k-1)-0.5)/real(nz-1)
  enddo

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
!    This sets the initial thickness (in m) of the layers.  The      !
!  thicknesses are set to insure that: 1.  each layer is at least an !
!  Angstrom thick, and 2.  the interfaces are where they should be   !
!  based on the resting depths and interface height perturbations,   !
!  as long at this doesn't interfere with 1.                         !
    eta1D(nz+1) = -1.0*G%bathyT(i,j)
    do k=nz,1,-1
      eta1D(K) = e0(K)
      if (eta1D(K) < (eta1D(K+1) + GV%Angstrom_z)) then
        eta1D(K) = eta1D(K+1) + GV%Angstrom_z
        h(i,j,k) = GV%Angstrom_z
      else
        h(i,j,k) = eta1D(K) - eta1D(K+1)
      endif
    enddo
  enddo ; enddo

end subroutine DOME_initialize_thickness
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> This subroutine sets the inverse restoration time (Idamp), and     !
!! the values towards which the interface heights and an arbitrary    !
!! number of tracers should be restored within each sponge. The       !
!! interface height is always subject to damping, and must always be  !
!! the first registered field.                                        !
subroutine DOME_initialize_sponges(G, GV, tv, PF, CSp)
  type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV !< The ocean's vertical grid structure.
  type(thermo_var_ptrs), intent(in) :: tv   !< A structure containing pointers to any available
               !!                 thermodynamic fields, including potential temperature and
               !!                 salinity or mixed layer density. Absent fields have NULL ptrs.
  type(param_file_type), intent(in) :: PF   !< A structure indicating the open file to
                                            !! parse for model parameter values.
  type(sponge_CS),       pointer    :: CSp  !< A pointer that is set to point to the control
                                            !! structure for this module.

  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! A temporary array for eta.
  real :: temp(SZI_(G),SZJ_(G),SZK_(G))  ! A temporary array for other variables. !
  real :: Idamp(SZI_(G),SZJ_(G))    ! The inverse damping rate, in s-1.

  real :: H0(SZK_(G))
  real :: min_depth
  real :: damp, e_dense, damp_new
  character(len=40)  :: mod = "DOME_initialize_sponges" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  eta(:,:,:) = 0.0 ; temp(:,:,:) = 0.0 ; Idamp(:,:) = 0.0

!  Here the inverse damping time, in s-1, is set. Set Idamp to 0     !
!  wherever there is no sponge, and the subroutines that are called  !
!  will automatically set up the sponges only where Idamp is positive!
!  and mask2dT is 1.                                                   !

!   Set up sponges for DOME configuration
  call get_param(PF, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)

  H0(1) = 0.0
  do k=2,nz ; H0(k) = -(real(k-1)-0.5)*G%max_depth/real(nz-1) ; enddo
  do i=is,ie; do j=js,je
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

! At this point, the DOME configuration is done. The following are here as a
! template for other configurations.

!  The remaining calls to set_up_sponge_field can be in any order. !
  if ( associated(tv%T) ) then
    call MOM_error(FATAL,"DOME_initialize_sponges is not set up for use with"//&
                         " a temperatures defined.")
    ! This should use the target values of T in temp.
    call set_up_sponge_field(temp, tv%T, G, nz, CSp)
    ! This should use the target values of S in temp.
    call set_up_sponge_field(temp, tv%S, G, nz, CSp)
  endif

end subroutine DOME_initialize_sponges

!> This subroutine sets the properties of flow at open boundary conditions.
!! This particular example is for the DOME inflow describe in Legg et al. 2006.
subroutine DOME_set_OBC_data(OBC, tv, G, GV, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC !< This open boundary condition type specifies
                                                !! whether, where, and what open boundary
                                                !! conditions are used.
  type(thermo_var_ptrs),      intent(in) :: tv  !< A structure containing pointers to any
                              !! available thermodynamic fields, including potential
                              !! temperature and salinity or mixed layer density. Absent
                              !! fields have NULL ptrs.
  type(ocean_grid_type),      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the open file
                              !! to parse for model parameter values.
  type(tracer_registry_type), pointer    :: tr_Reg !< Tracer registry.

  real, pointer, dimension(:,:,:) :: &
    OBC_T_v => NULL(), &    ! specify the values of T and S that should come
    OBC_S_v => NULL()       ! boundary conditions, in C and psu.
  ! The following variables are used to set the target temperature and salinity.
  real :: T0(SZK_(G)), S0(SZK_(G))
  real :: pres(SZK_(G))      ! An array of the reference pressure in Pa.
  real :: drho_dT(SZK_(G))   ! Derivative of density with temperature in kg m-3 K-1.                              !
  real :: drho_dS(SZK_(G))   ! Derivative of density with salinity in kg m-3 PSU-1.                             !
  real :: rho_guess(SZK_(G)) ! Potential density at T0 & S0 in kg m-3.
  ! The following variables are used to set up the transport in the DOME example.
  real :: tr_0, y1, y2, tr_k, rst, rsb, rc, v_k, lon_im1
  real :: D_edge            ! The thickness in m of the dense fluid at the
                            ! inner edge of the inflow.
  real :: g_prime_tot       ! The reduced gravity across all layers, m s-2.
  real :: Def_Rad           ! The deformation radius, based on fluid of
                            ! thickness D_edge, in the same units as lat.
  real :: Ri_trans          ! The shear Richardson number in the transition
                            ! region of the specified shear profile.
  character(len=40)  :: mod = "DOME_set_OBC_data" ! This subroutine's name.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! The following variables should be transformed into runtime parameters.
  D_edge = 300.0  ! The thickness of dense fluid in the inflow.
  Ri_trans = 1.0/3.0 ! The shear Richardson number in the transition region
                     ! region of the specified shear profile.

  if (.not.associated(OBC)) return

  g_prime_tot = (GV%g_Earth/GV%Rho0)*2.0
  Def_Rad = sqrt(D_edge*g_prime_tot) / (1.0e-4*1000.0)
  tr_0 = (-D_edge*sqrt(D_edge*g_prime_tot)*0.5e3*Def_Rad) * GV%m_to_H

  if (OBC%number_of_segments .ne. 1) then
    print *, 'Error in DOME OBC segment setup'
    return   !!! Need a better error message here
  endif
  segment => OBC%segment(1)
  if (.not. segment%on_pe) return

  do k=1,nz
    rst = -1.0
    if (k>1) rst = -1.0 + (real(k-1)-0.5)/real(nz-1)

    rsb = 0.0
    if (k<nz) rsb = -1.0 + (real(k-1)+0.5)/real(nz-1)
    rc = -1.0 + real(k-1)/real(nz-1)

    ! These come from assuming geostrophy and a constant Ri profile.
    y1 = (2.0*Ri_trans*rst + Ri_trans + 2.0)/(2.0 - Ri_trans)
    y2 = (2.0*Ri_trans*rsb + Ri_trans + 2.0)/(2.0 - Ri_trans)
    tr_k = tr_0 * (2.0/(Ri_trans*(2.0-Ri_trans))) * &
           ((log(y1)+1.0)/y1 - (log(y2)+1.0)/y2)
    v_k = -sqrt(D_edge*g_prime_tot)*log((2.0 + Ri_trans*(1.0 + 2.0*rc)) / &
                                        (2.0 - Ri_trans))
    if (k == nz)  tr_k = tr_k + tr_0 * (2.0/(Ri_trans*(2.0+Ri_trans))) * &
                                       log((2.0+Ri_trans)/(2.0-Ri_trans))
    ! New way
    isd = segment%HI%isd ; ied = segment%HI%ied
    JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
    do J=JsdB,JedB ; do i=isd,ied
      lon_im1 = 2.0*G%geoLonCv(i,J) - G%geoLonBu(I,J)
      segment%normal_trans(i,J,k) = tr_k * (exp(-2.0*(lon_im1 - 1000.0)/Def_Rad) -&
                                  exp(-2.0*(G%geoLonBu(I,J) - 1000.0)/Def_Rad))
      segment%normal_vel(i,J,k) = v_k * exp(-2.0*(G%geoLonCv(i,J) - 1000.0)/Def_Rad)
    enddo ; enddo
  enddo

  !   The inflow values of temperature and salinity also need to be set here if
  ! these variables are used.  The following code is just a naive example.
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

    allocate(OBC_T_v(isd:ied,JsdB:JedB,nz))
    do k=1,nz ; do J=JsdB,JedB ; do i=isd,ied
      OBC_T_v(i,J,k) = T0(k)
    enddo ; enddo ; enddo
    call add_tracer_OBC_values("T", tr_Reg, OBC_in_v=OBC_T_v)
  endif

end subroutine DOME_set_OBC_data

!> \class DOME_initialization
!!
!! The module configures the model for the "DOME" experiment.
!! DOME = Dynamics of Overflows and Mixing Experiment
end module DOME_initialization
