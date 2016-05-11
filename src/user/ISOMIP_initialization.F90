module ISOMIP_initialization
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
!*  The module configures the ISOMIP test case                         *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**
use MOM_ALE_sponge, only : ALE_sponge_CS, set_up_ALE_sponge_field, initialize_ALE_sponge
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, create_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use MOM_io, only : write_field, slasher, vardesc
use MOM_variables, only : thermo_var_ptrs, ocean_OBC_type
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA
implicit none ; private

#include <MOM_memory.h>

! -----------------------------------------------------------------------------
! Private (module-wise) parameters
! -----------------------------------------------------------------------------

character(len=40) :: mod = "ISOMIP_initialization" ! This module's name.

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public ISOMIP_initialize_topography
public ISOMIP_initialize_thickness
public ISOMIP_initialize_temperature_salinity 
public ISOMIP_initialize_sponges

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!> Initialization of topography
subroutine ISOMIP_initialize_topography(D, G, param_file, max_depth)
  type(ocean_grid_type), intent(in)           :: G
  real, intent(out), dimension(SZI_(G),SZJ_(G)) :: D
  type(param_file_type), intent(in)           :: param_file
  real,                  intent(in)           :: max_depth
! Arguments: D          - the bottom depth in m. Intent out.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets up the ISOMIP topography
  real :: min_depth ! The minimum and maximum depths in m.

! The following variables are used to set up the bathymetry in the ISOMIP example.
! check this paper: http://www.geosci-model-dev-discuss.net/8/9859/2015/gmdd-8-9859-2015.pdf

  real :: bmax            ! max depth of bedrock topography
  real :: b0,b2,b4,b6     ! first, second, third and fourth bedrock topography coeff
  real :: xbar           ! characteristic along-flow lenght scale of the bedrock
  real :: dc              ! depth of the trough compared with side walls
  real :: fc              ! characteristic width of the side walls of the channel
  real :: wc              ! half-width of the trough
  real :: ly              ! domain width (across ice flow)
  real :: bx, by, xtil    ! dummy vatiables
 
! G%ieg and G%jeg are the last indices in the global domain
!  real, dimension (G%ieg) :: xtil    ! eq. 3
!  real, dimension (G%ieg) :: bx            ! eq. 2
!  real, dimension (G%ieg,G%jeg) :: dummy1            ! dummy array
!  real, dimension (G%jeg) :: by            ! eq. 4

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "ISOMIP_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
 

  call MOM_mesg("  ISOMIP_initialization.F90, ISOMIP_initialize_topography: setting topography", 5)

  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)

! The following variables should be transformed into runtime parameters?
  bmax=720.0; b0=-150.0; b2=-728.8; b4=343.91; b6=-50.57
  xbar=300.0E3; dc=500.0; fc=4.0E3; wc=24.0E3; ly=80.0E3
  bx = 0.0; by = 0.0; xtil = 0.0

  do j=js,je ; do i=is,ie
       xtil = G%geoLonT(i,j)*1.0e3/xbar 
       bx = b0+b2*xtil**2 + b4*xtil**4 + b6*xtil**6
       by = (dc/(1.+exp(-2.*(G%geoLatT(i,j)*1.0e3- ly/2. - wc)/fc))) + &
            (dc/(1.+exp(2.*(G%geoLatT(i,j)*1.0e3- ly/2. + wc)/fc)))  
! depth is positive
        D(i,j) = -max(bx+by,-bmax)

    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
  enddo ; enddo

end subroutine ISOMIP_initialize_topography
! -----------------------------------------------------------------------------

!> Initialization of thicknesses
subroutine ISOMIP_initialize_thickness ( h, G, GV, param_file, tv )
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  real, intent(out), dimension(SZI_(G),SZJ_(G), SZK_(G)) :: h
  type(param_file_type), intent(in) :: param_file
  type(thermo_var_ptrs), intent(in) :: tv

! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields, including eq. of state
  real :: e0(SZK_(G)+1)     ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  integer :: i, j, k, is, ie, js, je, nz, tmp1
  real    :: x
  real    :: delta_h, rho_range
  real    :: min_thickness, s_sur, s_bot, t_sur, t_bot, rho_sur, rho_bot
  character(len=40) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("MOM_initialization.F90, initialize_thickness_uniform: setting thickness")

  call get_param(param_file,mod,"MIN_THICKNESS",min_thickness,'Minimum layer thickness',units='m',default=1.e-3)
  call get_param(param_file,mod,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
            default=DEFAULT_COORDINATE_MODE)
 
  select case ( coordinateMode(verticalCoordinate) )
    
  case ( REGRIDDING_LAYER, REGRIDDING_RHO ) ! Initial thicknesses for isopycnal coordinates
    call get_param(param_file, mod, "T_SUR",t_sur,'Temperature at the surface (interface)', default=-1.9)
    call get_param(param_file, mod, "S_SUF", s_sur, 'Salinity at the surface (interface)',  default=33.8)
    call get_param(param_file, mod, "T_BOT", t_bot, 'Temperature at the bottom (interface)', default=-1.9)
    call get_param(param_file, mod, "S_BOT", s_bot,'Salinity at the bottom (interface)', default=34.55)

    ! Compute min/max density using T_SUR/S_SUR and T_BOT/S_BOT
    call calculate_density(t_sur,s_sur,0.0,rho_sur,tv%eqn_of_state)
    write (*,*)'Surface density is:', rho_sur
    call calculate_density(t_bot,s_bot,0.0,rho_bot,tv%eqn_of_state)
    write (*,*)'Bottom density is:', rho_bot
    rho_range = rho_bot - rho_sur
    write (*,*)'Density range is:', rho_range
 
    ! Construct notional interface positions
    e0(1) = 0.
    do K=2,nz
      e0(k) = -G%max_depth * ( 0.5 * ( GV%Rlay(k-1) + GV%Rlay(k) ) - rho_sur ) / rho_range
      e0(k) = min( 0., e0(k) ) ! Bound by surface
      e0(k) = max( -G%max_depth, e0(k) ) ! Bound by possible deepest point in model
    enddo
    e0(nz+1) = -G%max_depth
    
    ! Calculate thicknesses
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -1.0*G%bathyT(i,j)
      do k=nz,1,-1
        eta1D(k) = e0(k)
        if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_z)) then
          eta1D(k) = eta1D(k+1) + GV%Angstrom_z
          h(i,j,k) = GV%Angstrom_z
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
    enddo ; enddo

  case ( REGRIDDING_ZSTAR )                       ! Initial thicknesses for z coordinates
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -1.0*G%bathyT(i,j)
      do k=nz,1,-1
        eta1D(k) =  -G%max_depth * real(k-1) / real(nz)
        if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
          eta1D(k) = eta1D(k+1) + min_thickness
          h(i,j,k) = min_thickness
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
   enddo ; enddo

  case ( REGRIDDING_SIGMA )             ! Initial thicknesses for sigma coordinates
    do j=js,je ; do i=is,ie
      delta_h = G%bathyT(i,j) / dfloat(nz)
      h(i,j,:) = delta_h
    end do ; end do 

  case default
      call MOM_error(FATAL,"isomip_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

end subroutine ISOMIP_initialize_thickness

!> Initial values for temperature and salinity
subroutine ISOMIP_initialize_temperature_salinity ( T, S, h, G, GV, param_file, &
                                                    eqn_of_state)
  type(ocean_grid_type),               intent(in)  :: G !< Ocean grid structure
  type(verticalGrid_type),                   intent(in) :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T !< Potential temperature (degC)
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: S !< Salinity (ppt)
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h !< Layer thickness (m or Pa)
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  type(EOS_type),                            pointer     :: eqn_of_state !< Equation of state structure
  ! Local variables

  integer   :: i, j, k, is, ie, js, je, nz
  real      :: x, ds, dt, rho_sur, rho_bot
  real      :: xi0, xi1, dxi, r, S_sur, T_sur, S_bot, T_bot, S_range, T_range
!  real    :: T_ref, T_Light, T_Dense, S_ref, S_Light, S_Dense, a1, frac_dense, k_frac, res_rat, k_light
  real      :: z          ! vertical position in z space
  character(len=40) :: verticalCoordinate, density_profile
real :: rho_light, delta_rho, z_tmp, rho_tmp
  
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call get_param(param_file,mod,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
            default=DEFAULT_COORDINATE_MODE)
  call get_param(param_file, mod, "T_SUR",t_sur,'Temperature at the surface (interface)', default=-1.9,&
                 do_not_log=.true.)
  call get_param(param_file, mod, "S_SUF", s_sur, 'Salinity at the surface (interface)',  default=33.8,&
                 do_not_log=.true.)
  call get_param(param_file, mod, "T_BOT", t_bot, 'Temperature at the bottom (interface)', default=-1.9,&
                 do_not_log=.true.)
  call get_param(param_file, mod, "S_BOT", s_bot,'Salinity at the bottom (interface)', default=34.55,&
                 do_not_log=.true.)

  call calculate_density(t_sur,s_sur,0.0,rho_sur,eqn_of_state)
  write (*,*)'Density in the surface layer:', rho_sur
  call calculate_density(t_bot,s_bot,0.0,rho_bot,eqn_of_state)
  write (*,*)'Density in the bottom layer::', rho_bot

  call get_param(param_file, mod, "LIGHTEST_DENSITY", rho_light)
  call get_param(param_file, mod, "DENSITY_RANGE", delta_rho)

  S_range = s_sur - s_bot
  T_range = t_sur - t_bot
  write(*,*)'s_bot, s_sur, S_range, t_bot, t_sur, T_range', s_bot, s_sur, S_range, t_bot, t_sur, T_range

  select case ( coordinateMode(verticalCoordinate) )

    case (  REGRIDDING_SIGMA, REGRIDDING_ZSTAR, REGRIDDING_RHO )
      S_range = S_range / G%max_depth ! Convert S_range into dS/dz
      T_range = T_range / G%max_depth ! Convert T_range into dT/dz
      do j=js,je ; do i=is,ie
        xi0 = -G%bathyT(i,j);
        do k = nz,1,-1
          xi0 = xi0 + 0.5 * h(i,j,k) ! Depth in middle of layer
          S(i,j,k) = S_sur + S_range * xi0
          T(i,j,k) = T_sur + T_range * xi0
          xi0 = xi0 + 0.5 * h(i,j,k) ! Depth at top of layer
        enddo
      enddo ; enddo
i=G%iec; j=G%jec
do k = 1,nz
  call calculate_density(T(i,j,k),S(i,j,k),0.0,rho_tmp,eqn_of_state)
  write(0,*) 'k,h,T,S,rho,Rlay',k,h(i,j,k),T(i,j,k),S(i,j,k),rho_tmp,GV%Rlay(k)
enddo
    
   case default
      call MOM_error(FATAL,"isomip_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

end subroutine ISOMIP_initialize_temperature_salinity

!> Sets up the the inverse restoration time (Idamp), and   !
! the values towards which the interface heights and an arbitrary    !
! number of tracers should be restored within each sponge.
subroutine ISOMIP_initialize_sponges(G,GV, tv, PF, CSp)
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  type(thermo_var_ptrs), intent(in) :: tv
  type(param_file_type), intent(in) :: PF
  type(ALE_sponge_CS),   pointer    :: CSp

! Arguments: G - The ocean's grid structure.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in)      PF - A structure indicating the open file to parse for
!                 model parameter values.
!  (in/out)  CSp - A pointer that is set to point to the control structure
!                  for this module

!  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! A temporary array for eta.
  real :: T(SZI_(G),SZJ_(G),SZK_(G))  ! A temporary array for temp 
  real :: S(SZI_(G),SZJ_(G),SZK_(G))  ! A temporary array for salt
  real :: RHO(SZI_(G),SZJ_(G),SZK_(G))  ! A temporary array for RHO
  real :: h(SZI_(G),SZJ_(G),SZK_(G))  ! A temporary array for thickness
  real :: Idamp(SZI_(G),SZJ_(G))    ! The inverse damping rate, in s-1.
  real      :: S_ref, T_ref;        ! Reference salinity and temerature within
                                    ! surface layer
  real      :: S_range, T_range;    ! Range of salinities and temperatures over the
                                    ! vertical


  real :: e0(SZK_(G))               ! The resting interface heights, in m, usually !
                                    ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)          ! Interface height relative to the sea surface !
                                    ! positive upward, in m.
  real :: min_depth, dummy1, z, Bmax
  real :: damp, rho_dummy
  character(len=40)  :: mod = "ISOMIP_initialize_sponges" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  T(:,:,:) = 0.0 ; S(:,:,:) = 0.0 ; Idamp(:,:) = 0.0; RHO(:,:,:) = 0.0
  S_ref = 33.8; S_range = 0.90
  T_ref = -1.9; T_range = 2.9
  Bmax = 720.0

!   Set up sponges for ISOMIP configuration
  call get_param(PF, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)

! GMM: set thickness, I am using same config. as the initial thickness for now
   do k=1,nz
    e0(k) = -G%max_depth * real(k-1) / real(nz)
   enddo

   do j=js,je ; do i=is,ie
     eta1D(nz+1) = -1.0*G%bathyT(i,j)
        do k=nz,1,-1
          eta1D(k) = e0(k)
          if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_z)) then
            eta1D(k) = eta1D(k+1) + GV%Angstrom_z
            h(i,j,k) = GV%Angstrom_z
          else
            h(i,j,k) = eta1D(k) - eta1D(k+1)
          endif
         enddo
    enddo;  enddo

!  Here the inverse damping time, in s-1, is set. Set Idamp to 0     !
!  wherever there is no sponge, and the subroutines that are called  !
!  will automatically set up the sponges only where Idamp is positive!
!  and mask2dT is 1.  

   do i=is,ie; do j=js,je
      if (G%geoLonT(i,j) >= 790.0 .AND. G%geoLonT(i,j) <= 800.0) then

! 1 / day
        dummy1=(G%geoLonT(i,j)-790.0)/(800.0-790.0)
!        damp = 1.0/10.0 * max(0.0,dummy1)
        damp = 1.0/0.04166 * max(0.0,dummy1)  ! one hour

      else ; damp=0.0
      endif

! convert to 1 / seconds
      if (G%bathyT(i,j) > min_depth) then
          Idamp(i,j) = damp/86400.0
      else ; Idamp(i,j) = 0.0 ; endif
   

   enddo ; enddo

   if (associated(CSp)) then
    call MOM_error(WARNING, "ISOMIP_initialize_sponges called with an associated "// &
                            "control structure.")
    return
  endif

!  This call sets up the damping rates and interface heights.
!  This sets the inverse damping timescale fields in the sponges.    !

  call initialize_ALE_sponge(Idamp,h, nz, G, PF, CSp)

! setup temp and salt at the sponge zone
      do j=js,je ; do i=is,ie
        do k = nz,1,-1
!          z = (G%bathyT(i,j)/(nz-1))* (k -1)
          z = (G%max_depth/(nz-1))* (k -1)
          S(i,j,k) = S_REF + (S_RANGE*z/Bmax)
          T(i,j,k) = T_REF + (T_RANGE*z/Bmax)
!          call calculate_density(T(i,j,k),S(i,j,k),0.0,rho_dummy,tv%eqn_of_state)
!          RHO(i,j,k) = rho_dummy
        enddo
      enddo ; enddo

!   Now register all of the fields which are damped in the sponge.   !
! By default, momentum is advected vertically within the sponge, but !
! momentum is typically not damped within the sponge.                !

! At this point, the ISOMIP configuration is done. The following are here as a
! template for other configurations.

!  The remaining calls to set_up_sponge_field can be in any order. !
  if ( associated(tv%T) ) then
      call set_up_ALE_sponge_field(T,G,tv%T,CSp)
  endif

  if ( associated(tv%S) ) then
    call set_up_ALE_sponge_field(S,G,tv%S,CSp)
  endif


end subroutine ISOMIP_initialize_sponges
! -----------------------------------------------------------------------------

end module ISOMIP_initialization
