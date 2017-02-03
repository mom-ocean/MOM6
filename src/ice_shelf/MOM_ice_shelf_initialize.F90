module MOM_ice_shelf_initialize
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                        *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

use MOM_grid, only : ocean_grid_type
use MOM_file_parser, only : get_param, read_param, log_param, param_file_type
use MOM_io, only: read_data, file_exists, slasher
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use user_shelf_init, only: USER_init_ice_thickness

implicit none ; private

#include <MOM_memory.h>
#ifdef SYMMETRIC_LAND_ICE
#  define GRID_SYM_ .true.
#  define NIMEMQ_IS_ NIMEMQS_
#  define NJMEMQ_IS_ NJMEMQS_
#  define ISUMSTART_INT_ CS%grid%iscq+1
#  define JSUMSTART_INT_ CS%grid%jscq+1
#else
#  define GRID_SYM_ .false.
#  define NIMEMQ_IS_ NIMEMQ_
#  define NJMEMQ_IS_ NJMEMQ_
#  define ISUMSTART_INT_ CS%grid%iscq
#  define JSUMSTART_INT_ CS%grid%jscq
#endif


!MJHpublic initialize_ice_shelf_boundary, initialize_ice_thickness
public initialize_ice_thickness

contains

subroutine initialize_ice_thickness (h_shelf, area_shelf_h, hmask, G, PF)

  real, intent(inout), dimension(:,:)          :: h_shelf, area_shelf_h, hmask
  type(ocean_grid_type), intent(in)            :: G
  type(param_file_type), intent(in)            :: PF

  character(len=40)  :: mod = "initialize_ice_thickness" ! This subroutine's name.
  character(len=200) :: config

  call get_param(PF, mod, "ICE_PROFILE_CONFIG", config, &
                 "This specifies how the initial ice profile is specified. \n"//&
                 "Valid values are: CHANNEL, FILE, and USER.", &
                 fail_if_missing=.true.)

  select case ( trim(config) )
  case ("CHANNEL"); call initialize_ice_thickness_channel (h_shelf, area_shelf_h, hmask, G, PF)
  case ("FILE");  call initialize_ice_thickness_from_file (h_shelf, area_shelf_h, hmask, G, PF)
  case ("USER");  call USER_init_ice_thickness (h_shelf, area_shelf_h, hmask, G, PF)
  case default ;  call MOM_error(FATAL,"MOM_initialize: "// &
    "Unrecognized ice profile setup "//trim(config))
  end select

end subroutine initialize_ice_thickness


subroutine initialize_ice_thickness_from_file (h_shelf, area_shelf_h, hmask, G, PF)

  real, intent(inout), dimension(:,:)          :: h_shelf, area_shelf_h, hmask
  type(ocean_grid_type), intent(in)            :: G
  type(param_file_type), intent(in)            :: PF

  !  This subroutine reads ice thickness and area from a file and puts it into
  !  h_shelf and area_shelf_h in m (and dimensionless) and updates hmask
  character(len=200) :: filename,thickness_file,inputdir ! Strings for file/path
  character(len=200) :: thickness_varname, area_varname! Variable name in file
  character(len=40)  :: mod = "initialize_ice_thickness_from_file" ! This subroutine's name.
  integer :: i, j, isc, jsc, iec, jec
  real :: len_sidestress, mask, udh

  call MOM_mesg("  MOM_ice_shelf_init_profile.F90, initialize_thickness_from_file: reading thickness")

  call get_param(PF, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(PF, mod, "ICE_THICKNESS_FILE", thickness_file, &
                 "The file from which the bathymetry is read.", &
                 default="ice_shelf_h.nc")
  call get_param(PF, mod, "LEN_SIDE_STRESS", len_sidestress, &
                 "position past which shelf sides are stress free.", &
                 default=0.0, units="axis_units")

  filename = trim(inputdir)//trim(thickness_file)
  call log_param(PF, mod, "INPUTDIR/THICKNESS_FILE", filename)
  call get_param(PF, mod, "ICE_THICKNESS_VARNAME", thickness_varname, &
                 "The name of the thickness variable in ICE_THICKNESS_FILE.", &
                 default="h_shelf")
  call get_param(PF, mod, "ICE_AREA_VARNAME", area_varname, &
                 "The name of the area variable in ICE_THICKNESS_FILE.", &
                 default="area_shelf_h")

  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topography_from_file: Unable to open "//trim(filename))

  call read_data(filename,trim(thickness_varname),h_shelf,domain=G%Domain%mpp_domain)
  call read_data(filename,trim(area_varname),area_shelf_h,domain=G%Domain%mpp_domain)

!  call get_param(PF, mod, "ICE_BOUNDARY_CONFIG", config, &
!                 "This specifies how the ice domain boundary is specified", &
!                 fail_if_missing=.true.)

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  do j=jsc,jec
    do i=isc,iec

      ! taper ice shelf in area where there is no sidestress -
      ! but do not interfere with hmask

      if ((G%geoLonCv(i,j) .gt. len_sidestress).and. &
          (len_sidestress .gt. 0.)) then
        udh = exp (-(G%geoLonCv(i,j)-len_sidestress)/5.0) * h_shelf(i,j)
        if (udh .le. 25.0) then
          h_shelf(i,j) = 0.0
          area_shelf_h (i,j) = 0.0
        else
          h_shelf(i,j) = udh
        endif
      endif

      ! update thickness mask

      if (area_shelf_h (i,j) .ge. G%areaT(i,j)) then
        hmask(i,j) = 1.
      elseif (area_shelf_h (i,j) .eq. 0.0) then
        hmask(i,j) = 0.
      elseif ((area_shelf_h(i,j) .gt. 0) .and. (area_shelf_h(i,j) .le. G%areaT(i,j))) then
        hmask(i,j) = 2.
      else
        call MOM_error(FATAL,mod// " AREA IN CELL OUT OF RANGE")
      endif
    enddo
  enddo


end subroutine initialize_ice_thickness_from_file


subroutine initialize_ice_thickness_channel (h_shelf, area_shelf_h, hmask, G, PF)

  real, intent(inout), dimension(:,:)          :: h_shelf, area_shelf_h, hmask
  type(ocean_grid_type), intent(in)            :: G
  type(param_file_type), intent(in)            :: PF

  character(len=40)  :: mod = "initialize_ice_shelf_thickness_channel" ! This subroutine's name.
  real :: max_draft, min_draft, flat_shelf_width, c1, slope_pos
  real :: edge_pos, shelf_slope_scale, Rho_ocean
  integer :: i, j, jsc, jec, jsd, jed, jedg, nyh, isc, iec, isd, ied
  integer :: j_off

  jsc = G%jsc ; jec = G%jec ; isc = G%isc ; iec = G%iec
  jsd = G%jsd ; jed = G%jed ; isd = G%isd ; ied = G%ied
  nyh = G%domain%njhalo ; jedg = G%domain%njglobal+nyh
  j_off = G%jdg_offset

  call MOM_mesg(mod//": setting thickness")

  call get_param(PF, mod, "SHELF_MAX_DRAFT", max_draft, &
                 units="m", default=1.0)
  call get_param(PF, mod, "SHELF_MIN_DRAFT", min_draft, &
                 units="m", default=1.0)
  call get_param(PF, mod, "FLAT_SHELF_WIDTH", flat_shelf_width, &
                 units="axis_units", default=0.0)
  call get_param(PF, mod, "SHELF_SLOPE_SCALE", shelf_slope_scale, &
                 units="axis_units", default=0.0)
  call get_param(PF, mod, "SHELF_EDGE_POS_0", edge_pos, &
                 units="axis_units", default=0.0)

  slope_pos = edge_pos - flat_shelf_width
  c1 = 0.0 ; if (shelf_slope_scale > 0.0) c1 = 1.0 / shelf_slope_scale

  do j=G%jsd,G%jed

  if (((j+j_off) <= jedg) .AND. ((j+j_off) >= nyh+1)) then

    do i=G%isc,G%iec

      if ((j.ge.jsc) .and. (j.le.jec)) then

        if (G%geoLonCu(i-1,j) >= edge_pos) then
        ! Everything past the edge is open ocean.
!            mass_shelf(i,j) = 0.0
          area_shelf_h(i,j) = 0.0
          hmask (i,j) = 0.0
          h_shelf (i,j) = 0.0
        else
          if (G%geoLonCu(i,j) > edge_pos) then
            area_shelf_h(i,j) = G%areaT(i,j) * (edge_pos - G%geoLonCu(i-1,j)) / &
                            (G%geoLonCu(i,j) - G%geoLonCu(i-1,j))
            hmask (i,j) = 2.0
          else
            area_shelf_h(i,j) = G%areaT(i,j)
            hmask (i,j) = 1.0
          endif

          if (G%geoLonT(i,j) > slope_pos) then
            h_shelf (i,j) = min_draft
!              mass_shelf(i,j) = Rho_ocean * min_draft
          else
!              mass_shelf(i,j) = Rho_ocean * (min_draft + &
!                 (CS%max_draft - CS%min_draft) * &
!                 min(1.0, (c1*(slope_pos - G%geoLonT(i,j)))**2) )
            h_shelf(i,j) = (min_draft + &
               (max_draft - min_draft) * &
               min(1.0, (c1*(slope_pos - G%geoLonT(i,j)))**2) )
          endif

        endif
      endif

      if ((i+G%idg_offset) .eq. G%domain%nihalo+1) then
        hmask(i-1,j) = 3.0
      endif

    enddo
  endif ; enddo

end subroutine initialize_ice_thickness_channel

!BEGIN MJH subroutine initialize_ice_shelf_boundary ( &
!      u_face_mask_boundary, &
!      v_face_mask_boundary, &
!      u_flux_boundary_values, &
!      v_flux_boundary_values, &
!      u_boundary_values, &
!      v_boundary_values, &
!      h_boundary_values, &
!      hmask, G, PF)

!   type(ocean_grid_type), intent(in)                 :: G
!   real, intent(inout), dimension(SZIB_(G),SZJ_(G))  :: u_face_mask_boundary, u_flux_boundary_values
!   real, intent(inout), dimension(SZI_(G),SZJB_(G))  :: v_face_mask_boundary, v_flux_boundary_values
!   real, intent(inout), dimension(SZIB_(G),SZJB_(G)) :: u_boundary_values, v_boundary_values
!   real, intent(inout), dimension(:,:)               :: hmask, h_boundary_values
!   type(param_file_type), intent(in)                 :: PF

!   character(len=40)  :: mod = "initialize_ice_shelf_boundary" ! This subroutine's name.
!   character(len=200) :: config
!   logical flux_bdry

!   call get_param(PF, mod, "ICE_BOUNDARY_CONFIG", config, &
!                  "This specifies how the ice domain boundary is specified. \n"//&
!                  "valid values include CHANNEL, FILE and USER.", &
!                  fail_if_missing=.true.)
!   call get_param(PF, mod, "ICE_BOUNDARY_FLUX_CONDITION", flux_bdry, &
!                  "This specifies whether mass input is a dirichlet or \n"//&
!                  "flux condition", default=.true.)

!   select case ( trim(config) )
!     case ("CHANNEL");
!       call initialize_ice_shelf_boundary_channel(u_face_mask_boundary, &
!         v_face_mask_boundary, u_flux_boundary_values, v_flux_boundary_values, &
!         u_boundary_values, v_boundary_values, h_boundary_values, hmask, G, &
!         flux_bdry, PF)
!     case ("FILE");  call MOM_error(FATAL,"MOM_initialize: "// &
!       "Unrecognized topography setup "//trim(config))
!     case ("USER");  call MOM_error(FATAL,"MOM_initialize: "// &
!       "Unrecognized topography setup "//trim(config))
!     case default ;  call MOM_error(FATAL,"MOM_initialize: "// &
!       "Unrecognized topography setup "//trim(config))
!   end select

! end subroutine initialize_ice_shelf_boundary

! subroutine initialize_ice_shelf_boundary_channel ( &
!   u_face_mask_boundary, &
!   v_face_mask_boundary, &
!   u_flux_boundary_values, &
!   v_flux_boundary_values, &
!   u_boundary_values, &
!   v_boundary_values, &
!   h_boundary_values, &
!   hmask, &
!   G, flux_bdry, PF )

!   type(ocean_grid_type), intent(in)                 :: G
!   real, dimension(SZIB_(G),SZJ_(G)),  intent(inout) :: u_face_mask_boundary, u_flux_boundary_values
!   real, dimension(SZI_(G),SZJB_(G)),  intent(inout) :: v_face_mask_boundary, v_flux_boundary_values
!   real, dimension(SZIB_(G),SZJB_(G)), intent(inout) :: u_boundary_values, v_boundary_values
!   real, dimension(:,:), intent(inout)               :: h_boundary_values, hmask
!   logical, intent(in)                               :: flux_bdry
!   type (param_file_type), intent(in)                :: PF

!   character(len=40)  :: mod = "initialize_ice_shelf_boundary_channel" ! This subroutine's name.
!   integer :: i, j, isd, jsd, is, js, iegq, jegq, giec, gjec, gisc, gjsc, isc, jsc, iec, jec, ied, jed
!   real                                                  :: lenlat, input_thick, input_flux, len_stress

!   call get_param(PF, mod, "LENLAT", lenlat, fail_if_missing=.true.)

!   call get_param(PF, mod, "INPUT_FLUX_ICE_SHELF", input_flux, &
!                  "volume flux at upstream boundary", &
!                  units="m2 s-1", default=0.)
!   call get_param(PF, mod, "INPUT_THICK_ICE_SHELF", input_thick, &
!                  "flux thickness at upstream boundary", &
!                  units="m", default=1000.)
!   call get_param(PF, mod, "LEN_SIDE_STRESS", len_stress, &
!                  "maximum position of no-flow condition in along-flow direction", &
!                  units="km", default=0.)

!   call MOM_mesg(mod//": setting boundary")

!   isd = G%isd ; ied = G%ied
!   jsd = G%jsd ; jed = G%jed
!   isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
!   gisc = G%Domain%nihalo ; gjsc = G%Domain%njhalo
!   giec = G%Domain%niglobal+gisc ; gjec = G%Domain%njglobal+gjsc

!   do j=jsd,jed
!     do i=isd,ied

!       ! upstream boundary - set either dirichlet or flux condition

!       if ((i+G%idg_offset) .eq. G%domain%nihalo+1) then
!         if (flux_bdry) then
!           u_face_mask_boundary (i-1,j) = 4.0
!           u_flux_boundary_values (i-1,j) = input_flux
!         else
!           hmask(i-1,j) = 3.0
!           h_boundary_values (i-1,j) = input_thick
!           u_face_mask_boundary (i-1,j) = 3.0
!           u_boundary_values (i-1,j-1) = (1 - ((G%geoLatBu(i-1,j-1) - 0.5*lenlat)*2./lenlat)**2) * &
!                   1.5 * input_flux / input_thick
!           u_boundary_values (i-1,j) = (1 - ((G%geoLatBu(i-1,j) - 0.5*lenlat)*2./lenlat)**2) * &
!                   1.5 * input_flux / input_thick
!         endif
!       endif

!       ! side boundaries: no flow

!       if (G%jdg_offset+j .eq. gjsc+1) then !bot boundary
!         if (len_stress .eq. 0. .OR. G%geoLonCv(i,j-1) .le. len_stress) then
!           v_face_mask_boundary (i,j-1) = 0.
!         else
!           v_face_mask_boundary (i,j-1) = 1.
!         endif
!       elseif (G%jdg_offset+j .eq. gjec) then !top boundary
!         if (len_stress .eq. 0. .OR. G%geoLonCv(i,j-1) .le. len_stress) then
!           v_face_mask_boundary (i,j) = 0.
!         else
!           v_face_mask_boundary (i,j) = 1.
!         endif
!       endif

!       ! downstream boundary - CFBC

!       if (i+G%idg_offset .eq. giec) then
!         u_face_mask_boundary(i,j) = 2.0
!       endif

!     enddo
!   enddo

!END MJH end subroutine initialize_ice_shelf_boundary_channel

end module MOM_ice_shelf_initialize
