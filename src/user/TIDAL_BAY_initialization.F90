module TIDAL_BAY_initialization
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

use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_open_boundary, only : open_boundary_query, set_Flather_positions
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public TIDAL_BAY_set_OBC_positions
public TIDAL_BAY_set_OBC_data

contains

!> Set the positions of the open boundary needed for the TIDAL_BAY experiment.
subroutine TIDAL_BAY_set_OBC_positions(G, param_file, OBC)
  type(dyn_horgrid_type),     intent(inout) :: G   !< Grid structure.
  type(param_file_type),      intent(in)    :: param_file !< Parameter file handle.
  type(ocean_OBC_type),       pointer       :: OBC !< Open boundary control structure.
  ! Local variables
  character(len=40)  :: mod = "TIDAL_BAY_set_OBC_positions" ! This subroutine's name.
  integer :: i, j

  if (.not.associated(OBC)) call MOM_error(FATAL, &
           "TIDAL_BAY_initialization, TIDAL_BAY_set_OBC_positions: OBC type was not allocated!")

  ! This isn't called when APPLY_OBC_U is requested.
  if (open_boundary_query(OBC, apply_orig_Flather=.true.)) then
    call set_Flather_positions(G, OBC)
  endif
  if (OBC%apply_OBC_u) then
    ! Set where u points are determined by OBCs.
    allocate(OBC%OBC_mask_u(G%isd:G%ied,G%JsdB:G%JedB)) ; OBC%OBC_mask_u(:,:) = .false.
    do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
      if ((G%geoLonCv(i,J) > 1000.0) .and. (G%geoLonCv(i,J)  < 1100.0) .and. &
          (abs(G%geoLatCv(i,J) - G%gridLatB(G%JegB)) < 0.1)) then
        OBC%OBC_mask_u(i,J) = .true.
      endif
    enddo ; enddo
  endif
  if (OBC%apply_OBC_v) then
    ! Set where u points are determined by OBCs.
    !allocate(OBC_mask_u(IsdB:IedB,jsd:jed)) ; OBC_mask_u(:,:) = .false.
    call MOM_error(FATAL,"TIDAL_BAY_initialization, TIDAL_BAY_set_OBC_positions: "//&
                   "APPLY_OBC_U=True is not coded for the TIDAL_BAY experiment")
  endif
end subroutine TIDAL_BAY_set_OBC_positions

!> This subroutine sets the properties of flow at open boundary conditions.
!! This particular example is for the TIDAL_BAY inflow describe in Legg et al. 2006.
subroutine TIDAL_BAY_set_OBC_data(OBC, G, GV, param_file)
  type(ocean_OBC_type),       pointer    :: OBC !< This open boundary condition type specifies
                                                !! whether, where, and what open boundary
                                                !! conditions are used.
  type(ocean_grid_type),      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the open file
                              !! to parse for model parameter values.

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
  ! The following variables are used to set up the transport in the TIDAL_BAY example.
  real :: tr_0, y1, y2, tr_k, rst, rsb, rc, v_k, lon_im1
  real :: D_edge            ! The thickness in m of the dense fluid at the
                            ! inner edge of the inflow.
  real :: g_prime_tot       ! The reduced gravity across all layers, m s-2.
  real :: Def_Rad           ! The deformation radius, based on fluid of
                            ! thickness D_edge, in the same units as lat.
  real :: Ri_trans          ! The shear Richardson number in the transition
                            ! region of the specified shear profile.
  character(len=40)  :: mod = "TIDAL_BAY_set_OBC_data" ! This subroutine's name.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! The following variables should be transformed into runtime parameters.
  D_edge = 300.0  ! The thickness of dense fluid in the inflow.
  Ri_trans = 1.0/3.0 ! The shear Richardson number in the transition region
                     ! region of the specified shear profile.

  if (.not.associated(OBC)) return
  if (.not.(OBC%apply_OBC_u .or. OBC%apply_OBC_v)) return

  if (OBC%apply_OBC_u) then
    allocate(OBC%u(IsdB:IedB,jsd:jed,nz)) ; OBC%u(:,:,:) = 0.0
    allocate(OBC%uh(IsdB:IedB,jsd:jed,nz)) ; OBC%uh(:,:,:) = 0.0
    allocate(OBC%OBC_kind_u(IsdB:IedB,jsd:jed)) ; OBC%OBC_kind_u(:,:) = OBC_NONE
    allocate(OBC%OBC_direction_u(IsdB:IedB,jsd:jed)) ; OBC%OBC_direction_u(:,:) = OBC_NONE
    do j=jsd,jed ; do I=IsdB,IedB
      if (OBC%OBC_mask_u(I,j)) OBC%OBC_kind_u(I,j) = OBC_SIMPLE
    enddo ; enddo
  endif
  if (OBC%apply_OBC_v) then
    allocate(OBC%v(isd:ied,JsdB:JedB,nz)) ; OBC%v(:,:,:) = 0.0
    allocate(OBC%vh(isd:ied,JsdB:JedB,nz)) ; OBC%vh(:,:,:) = 0.0
    allocate(OBC%OBC_kind_v(isd:ied,JsdB:JedB)) ; OBC%OBC_kind_v(:,:) = OBC_NONE
    allocate(OBC%OBC_direction_v(isd:ied,JsdB:JedB)) ; OBC%OBC_direction_v(:,:) = OBC_NONE
    do J=JsdB,JedB ; do i=isd,ied
      if (OBC%OBC_mask_v(i,J)) OBC%OBC_kind_v(i,J) = OBC_SIMPLE
    enddo ; enddo
  endif

  if (OBC%apply_OBC_v) then
    g_prime_tot = (GV%g_Earth/GV%Rho0)*2.0
    Def_Rad = sqrt(D_edge*g_prime_tot) / (1.0e-4*1000.0)
    tr_0 = (-D_edge*sqrt(D_edge*g_prime_tot)*0.5e3*Def_Rad) * GV%m_to_H

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
      do J=JsdB,JedB ; do i=isd,ied
        if (OBC%OBC_mask_v(i,J)) then
          ! This needs to be unneccesarily complicated without symmetric memory.
          lon_im1 = 2.0*G%geoLonCv(i,J) - G%geoLonBu(I,J)
          ! if (isd > IsdB) lon_im1 = G%geoLonBu(I-1,J)
          OBC%vh(i,J,k) = tr_k * (exp(-2.0*(lon_im1 - 1000.0)/Def_Rad) -&
                                exp(-2.0*(G%geoLonBu(I,J) - 1000.0)/Def_Rad))
          OBC%v(i,J,k) = v_k * exp(-2.0*(G%geoLonCv(i,J) - 1000.0)/Def_Rad)
        else
          OBC%vh(i,J,k) = 0.0 ; OBC%v(i,J,k) = 0.0
        endif
      enddo ; enddo
    enddo
  endif

  if (OBC%apply_OBC_u) then
    do k=1,nz ; do j=jsd,jed ; do I=IsdB,IedB
      if (OBC%OBC_mask_u(I,j)) then
        ! An appropriate expression for the zonal inflow velocities and
        ! transports should go here.
        OBC%uh(I,j,k) = 0.0 * GV%m_to_H ; OBC%u(I,j,k) = 0.0
      else
        OBC%uh(I,j,k) = 0.0 ; OBC%u(I,j,k) = 0.0
      endif
    enddo ; enddo ; enddo
  endif

  !   The inflow values of temperature and salinity also need to be set here if
  ! these variables are used.  The following code is just a naive example.
  if (OBC%apply_OBC_u .or. OBC%apply_OBC_v) then
  endif

end subroutine TIDAL_BAY_set_OBC_data

!> \class TIDAL_BAY_initialization
!!
!! The module configures the model for the "TIDAL_BAY" experiment.
!! TIDAL_BAY = Dynamics of Overflows and Mixing Experiment
end module TIDAL_BAY_initialization
