!> Initialize ice shelf variables
module MOM_ice_shelf_initialize

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_grid, only : ocean_grid_type
use MOM_array_transform,      only : rotate_array
use MOM_hor_index,  only : hor_index_type
use MOM_file_parser, only : get_param, read_param, log_param, param_file_type
use MOM_io, only: MOM_read_data, file_exists, slasher
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_unit_scaling, only : unit_scale_type
use user_shelf_init, only: USER_init_ice_thickness

implicit none ; private

#include <MOM_memory.h>

!MJHpublic initialize_ice_shelf_boundary, initialize_ice_thickness
public initialize_ice_thickness

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> Initialize ice shelf thickness
subroutine initialize_ice_thickness(h_shelf, area_shelf_h, hmask, G, G_in, US, PF, rotate_index, turns)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(ocean_grid_type), intent(in)    :: G_in    !< The ocean's unrotated grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters
  logical, intent(in), optional        :: rotate_index !< If true, this is a rotation test
  integer, intent(in), optional        :: turns !< Number of turns for rotation test

  integer :: i, j
  character(len=40)  :: mdl = "initialize_ice_thickness" ! This subroutine's name.
  character(len=200) :: config
  logical :: rotate = .false.
  real, allocatable, dimension(:,:) :: tmp1_2d ! Temporary array for storing ice shelf input data
  real, allocatable, dimension(:,:) :: tmp2_2d ! Temporary array for storing ice shelf input data
  real, allocatable, dimension(:,:) :: tmp3_2d ! Temporary array for storing ice shelf input data

  call get_param(PF, mdl, "ICE_PROFILE_CONFIG", config, &
                 "This specifies how the initial ice profile is specified. "//&
                 "Valid values are: CHANNEL, FILE, and USER.", &
                 fail_if_missing=.true.)

  if (PRESENT(rotate_index)) rotate=rotate_index

  if (rotate) then
     allocate(tmp1_2d(G_in%isd:G_in%ied,G_in%jsd:G_in%jed));tmp1_2d(:,:)=0.0
     allocate(tmp2_2d(G_in%isd:G_in%ied,G_in%jsd:G_in%jed));tmp2_2d(:,:)=0.0
     allocate(tmp3_2d(G_in%isd:G_in%ied,G_in%jsd:G_in%jed));tmp3_2d(:,:)=0.0
     select case ( trim(config) )
     case ("CHANNEL"); call initialize_ice_thickness_channel (tmp1_2d, tmp2_2d, tmp3_2d, G_in, US, PF)
     case ("FILE");  call initialize_ice_thickness_from_file (tmp1_2d, tmp2_2d, tmp3_2d, G_in, US, PF)
     case ("USER");  call USER_init_ice_thickness (tmp1_2d, tmp2_2d, tmp3_2d, G_in, US, PF)
     case default ;  call MOM_error(FATAL,"MOM_initialize: "// &
          "Unrecognized ice profile setup "//trim(config))
     end select
     call rotate_array(tmp1_2d,turns, h_shelf)
     call rotate_array(tmp2_2d,turns, area_shelf_h)
     call rotate_array(tmp3_2d,turns, hmask)
     deallocate(tmp1_2d,tmp2_2d,tmp3_2d)
  else
     select case ( trim(config) )
     case ("CHANNEL"); call initialize_ice_thickness_channel (h_shelf, area_shelf_h, hmask, G, US, PF)
     case ("FILE");  call initialize_ice_thickness_from_file (h_shelf, area_shelf_h, hmask, G, US, PF)
     case ("USER");  call USER_init_ice_thickness (h_shelf, area_shelf_h, hmask, G, US, PF)
     case default ;  call MOM_error(FATAL,"MOM_initialize: "// &
          "Unrecognized ice profile setup "//trim(config))
     end select
  endif

end subroutine initialize_ice_thickness

!> Initialize ice shelf thickness from file
subroutine initialize_ice_thickness_from_file(h_shelf, area_shelf_h, hmask, G, US, PF)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters

  !  This subroutine reads ice thickness and area from a file and puts it into
  !  h_shelf [Z ~> m] and area_shelf_h [L2 ~> m2] (and dimensionless) and updates hmask
  character(len=200) :: filename,thickness_file,inputdir ! Strings for file/path
  character(len=200) :: thickness_varname, area_varname  ! Variable name in file
  character(len=40)  :: mdl = "initialize_ice_thickness_from_file" ! This subroutine's name.
  integer :: i, j, isc, jsc, iec, jec
  real :: len_sidestress, mask, udh

  call MOM_mesg("Initialize_ice_thickness_from_file: reading thickness")

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(PF, mdl, "ICE_THICKNESS_FILE", thickness_file, &
                 "The file from which the bathymetry is read.", &
                 default="ice_shelf_h.nc")
  call get_param(PF, mdl, "LEN_SIDE_STRESS", len_sidestress, &
                 "position past which shelf sides are stress free.", &
                 default=0.0, units="axis_units")

  filename = trim(inputdir)//trim(thickness_file)
  call log_param(PF, mdl, "INPUTDIR/THICKNESS_FILE", filename)
  call get_param(PF, mdl, "ICE_THICKNESS_VARNAME", thickness_varname, &
                 "The name of the thickness variable in ICE_THICKNESS_FILE.", &
                 default="h_shelf")
  call get_param(PF, mdl, "ICE_AREA_VARNAME", area_varname, &
                 "The name of the area variable in ICE_THICKNESS_FILE.", &
                 default="area_shelf_h")

  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topography_from_file: Unable to open "//trim(filename))
  call MOM_read_data(filename, trim(thickness_varname), h_shelf, G%Domain, scale=US%m_to_Z)
  call MOM_read_data(filename,trim(area_varname), area_shelf_h, G%Domain, scale=US%m_to_L**2)

!  call get_param(PF, mdl, "ICE_BOUNDARY_CONFIG", config, &
!                 "This specifies how the ice domain boundary is specified", &
!                 fail_if_missing=.true.)

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  do j=jsc,jec
    do i=isc,iec

      ! taper ice shelf in area where there is no sidestress -
      ! but do not interfere with hmask

      if ((G%geoLonCv(i,j) > len_sidestress).and. &
          (len_sidestress > 0.)) then
        udh = exp(-(G%geoLonCv(i,j)-len_sidestress)/5.0) * h_shelf(i,j)
        if (udh <= 25.0) then
          h_shelf(i,j) = 0.0
          area_shelf_h(i,j) = 0.0
        else
          h_shelf(i,j) = udh
        endif
      endif

      ! update thickness mask

      if (area_shelf_h (i,j) >= G%areaT(i,j)) then
        hmask(i,j) = 1.
      elseif (area_shelf_h (i,j) == 0.0) then
        hmask(i,j) = 0.
      elseif ((area_shelf_h(i,j) > 0) .and. (area_shelf_h(i,j) <= G%areaT(i,j))) then
        hmask(i,j) = 2.
      else
        call MOM_error(FATAL,mdl// " AREA IN CELL OUT OF RANGE")
      endif
    enddo
  enddo

end subroutine initialize_ice_thickness_from_file

!> Initialize ice shelf thickness for a channel configuration
subroutine initialize_ice_thickness_channel(h_shelf, area_shelf_h, hmask, G, US, PF)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters

  character(len=40)  :: mdl = "initialize_ice_shelf_thickness_channel" ! This subroutine's name.
  real :: max_draft, min_draft, flat_shelf_width, c1, slope_pos
  real :: edge_pos, shelf_slope_scale, Rho_ocean
  integer :: i, j, jsc, jec, jsd, jed, jedg, nyh, isc, iec, isd, ied
  integer :: j_off

  jsc = G%jsc ; jec = G%jec ; isc = G%isc ; iec = G%iec
  jsd = G%jsd ; jed = G%jed ; isd = G%isd ; ied = G%ied
  nyh = G%domain%njhalo ; jedg = G%domain%njglobal+nyh
  j_off = G%jdg_offset

  call MOM_mesg(mdl//": setting thickness")

  call get_param(PF, mdl, "SHELF_MAX_DRAFT", max_draft, &
                 units="m", default=1.0, scale=US%m_to_Z)
  call get_param(PF, mdl, "SHELF_MIN_DRAFT", min_draft, &
                 units="m", default=1.0, scale=US%m_to_Z)
  call get_param(PF, mdl, "FLAT_SHELF_WIDTH", flat_shelf_width, &
                 units="axis_units", default=0.0)
  call get_param(PF, mdl, "SHELF_SLOPE_SCALE", shelf_slope_scale, &
                 units="axis_units", default=0.0)
  call get_param(PF, mdl, "SHELF_EDGE_POS_0", edge_pos, &
                 units="axis_units", default=0.0)
!  call get_param(param_file, mdl, "RHO_0", Rho_ocean, &
!                 "The mean ocean density used with BOUSSINESQ true to "//&
!                 "calculate accelerations and the mass for conservation "//&
!                 "properties, or with BOUSSINSEQ false to convert some "//&
!                 "parameters from vertical units of m to kg m-2.", &
!                 units="kg m-3", default=1035.0, scale=US%Z_to_m)

  slope_pos = edge_pos - flat_shelf_width
  c1 = 0.0 ; if (shelf_slope_scale > 0.0) c1 = 1.0 / shelf_slope_scale


  do j=G%jsd,G%jed

  if (((j+j_off) <= jedg) .AND. ((j+j_off) >= nyh+1)) then

    do i=G%isc,G%iec

      if ((j >= jsc) .and. (j <= jec)) then

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
            h_shelf(i,j) = min_draft
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

      if ((i+G%idg_offset) == G%domain%nihalo+1) then
        hmask(i-1,j) = 3.0
      endif

    enddo
  endif ; enddo

end subroutine initialize_ice_thickness_channel

!BEGIN MJH
! subroutine initialize_ice_shelf_boundary(u_face_mask_bdry, v_face_mask_bdry, &
!                u_flux_bdry_val, v_flux_bdry_val, u_bdry_val, v_bdry_val, h_bdry_val, &
!                hmask, G, US, PF )

!   type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
!   real, dimension(SZIB_(G),SZJ_(G)), &
!                          intent(inout) :: u_face_mask_bdry !< A boundary-type mask at C-grid u faces
!   real, dimension(SZIB_(G),SZJ_(G)), &
!                          intent(inout) :: u_flux_bdry_val  !< The boundary thickness flux through
!                                                      !! C-grid u faces [L Z T-1 ~> m2 s-1].
!   real, dimension(SZI_(G),SZJB_(G)), &
!                          intent(inout) :: v_face_mask_bdry !< A boundary-type mask at C-grid v faces
!   real, dimension(SZI_(G),SZJB_(G)), &
!                          intent(inout) :: v_flux_bdry_val  !< The boundary thickness flux through
!                                                      !! C-grid v faces [L Z T-1 ~> m2 s-1].
!   real, dimension(SZIB_(G),SZJB_(G)), &
!                          intent(inout) :: u_bdry_val !< The zonal ice shelf velocity at open
!                                                      !! boundary vertices [L T-1 ~> m s-1].
!   real, dimension(SZIB_(G),SZJB_(G)), &
!                          intent(inout) :: v_bdry_val !< The meridional ice shelf velocity at open
!                                                      !! boundary vertices [L T-1 ~> m s-1].
!   real, dimension(SZDI_(G),SZDJ_(G)), &
!                          intent(inout) :: h_bdry_val !< The ice shelf thickness at open boundaries [Z ~> m]
!   real, dimension(SZDI_(G),SZDJ_(G)), &
!                          intent(inout) :: hmask !< A mask indicating which tracer points are
!                                              !! partly or fully covered by an ice-shelf
!   type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
!   type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters

!   character(len=40)  :: mdl = "initialize_ice_shelf_boundary" ! This subroutine's name.
!   character(len=200) :: config
!   logical flux_bdry

!   call get_param(PF, mdl, "ICE_BOUNDARY_CONFIG", config, &
!                  "This specifies how the ice domain boundary is specified. "//&
!                  "valid values include CHANNEL, FILE and USER.", &
!                  fail_if_missing=.true.)
!   call get_param(PF, mdl, "ICE_BOUNDARY_FLUX_CONDITION", flux_bdry, &
!                  "This specifies whether mass input is a dirichlet or "//&
!                  "flux condition", default=.true.)

!   select case ( trim(config) )
!     case ("CHANNEL")
!       call initialize_ice_shelf_boundary_channel(u_face_mask_bdry, &
!         v_face_mask_bdry, u_flux_bdry_val, v_flux_bdry_val, &
!         u_bdry_val, v_bdry_val, h_bdry_val, hmask, G, &
!         flux_bdry, PF)
!     case ("FILE");  call MOM_error(FATAL,"MOM_initialize: "// &
!       "Unrecognized topography setup "//trim(config))
!     case ("USER");  call MOM_error(FATAL,"MOM_initialize: "// &
!       "Unrecognized topography setup "//trim(config))
!     case default ;  call MOM_error(FATAL,"MOM_initialize: "// &
!       "Unrecognized topography setup "//trim(config))
!   end select

! end subroutine initialize_ice_shelf_boundary

! subroutine initialize_ice_shelf_boundary_channel(u_face_mask_bdry, v_face_mask_bdry, &
!                u_flux_bdry_val, v_flux_bdry_val, u_bdry_val, v_bdry_val, h_bdry_val, &
!                hmask, G, flux_bdry, US, PF )

!   type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
!   real, dimension(SZIB_(G),SZJ_(G)), &
!                          intent(inout) :: u_face_mask_bdry !< A boundary-type mask at C-grid u faces
!   real, dimension(SZIB_(G),SZJ_(G)), &
!                          intent(inout) :: u_flux_bdry_val  !< The boundary thickness flux through
!                                                      !! C-grid u faces [L Z T-1 ~> m2 s-1].
!   real, dimension(SZI_(G),SZJB_(G)), &
!                          intent(inout) :: v_face_mask_bdry !< A boundary-type mask at C-grid v faces
!   real, dimension(SZI_(G),SZJB_(G)), &
!                          intent(inout) :: v_flux_bdry_val  !< The boundary thickness flux through
!                                                      !! C-grid v faces [L Z T-1 ~> m2 s-1].
!   real, dimension(SZIB_(G),SZJB_(G)), &
!                          intent(inout) :: u_bdry_val !< The zonal ice shelf velocity at open
                                                       !! boundary vertices [L T-1 ~> m s-1].
!   real, dimension(SZIB_(G),SZJB_(G)), &
!                          intent(inout) :: v_bdry_val !< The meridional ice shelf velocity at open
                                                       !! boundary vertices [L T-1 ~> m s-1].
!   real, dimension(SZDI_(G),SZDJ_(G)), &
!                          intent(inout) :: h_bdry_val !< The ice shelf thickness at open boundaries [Z ~> m]
!   real, dimension(SZDI_(G),SZDJ_(G)), &
!                          intent(inout) :: hmask !< A mask indicating which tracer points are
!                                              !! partly or fully covered by an ice-shelf
!   logical,               intent(in)    :: flux_bdry !< If true, use mass fluxes as the boundary value.
!   type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
!   type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters

!   character(len=40)  :: mdl = "initialize_ice_shelf_boundary_channel" ! This subroutine's name.
!   integer :: i, j, isd, jsd, is, js, iegq, jegq, giec, gjec, gisc, gjsc, isc, jsc, iec, jec, ied, jed
!   real    :: input_thick ! The input ice shelf thickness [Z ~> m]
!   real    :: input_flux  ! The input ice flux per unit length [L Z T-1 ~> m2 s-1]
!   real    :: lenlat, len_stress

!   call get_param(PF, mdl, "LENLAT", lenlat, fail_if_missing=.true.)

!   call get_param(PF, mdl, "INPUT_FLUX_ICE_SHELF", input_flux, &
!                  "volume flux at upstream boundary", &
!                  units="m2 s-1", default=0., scale=US%m_s_to_L_T*US%m_to_Z)
!   call get_param(PF, mdl, "INPUT_THICK_ICE_SHELF", input_thick, &
!                  "flux thickness at upstream boundary", &
!                  units="m", default=1000., scale=US%m_to_Z)
!   call get_param(PF, mdl, "LEN_SIDE_STRESS", len_stress, &
!                  "maximum position of no-flow condition in along-flow direction", &
!                  units="km", default=0.)

!   call MOM_mesg(mdl//": setting boundary")

!   isd = G%isd ; ied = G%ied
!   jsd = G%jsd ; jed = G%jed
!   isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
!   gisc = G%Domain%nihalo ; gjsc = G%Domain%njhalo
!   giec = G%Domain%niglobal+gisc ; gjec = G%Domain%njglobal+gjsc

!   do j=jsd,jed
!     do i=isd,ied

!       ! upstream boundary - set either dirichlet or flux condition

!       if ((i+G%idg_offset) == G%domain%nihalo+1) then
!         if (flux_bdry) then
!           u_face_mask_bdry(i-1,j) = 4.0
!           u_flux_bdry_val(i-1,j) = input_flux
!         else
!           hmask(i-1,j) = 3.0
!           h_bdry_val(i-1,j) = input_thick
!           u_face_mask_bdry(i-1,j) = 3.0
!           u_bdry_val(i-1,j-1) = (1 - ((G%geoLatBu(i-1,j-1) - 0.5*lenlat)*2./lenlat)**2) * &
!                   1.5 * input_flux / input_thick
!           u_bdry_val(i-1,j) = (1 - ((G%geoLatBu(i-1,j) - 0.5*lenlat)*2./lenlat)**2) * &
!                   1.5 * input_flux / input_thick
!         endif
!       endif

!       ! side boundaries: no flow

!       if (G%jdg_offset+j == gjsc+1) then !bot boundary
!         if (len_stress == 0. .OR. G%geoLonCv(i,j-1) <= len_stress) then
!           v_face_mask_bdry(i,j-1) = 0.
!         else
!           v_face_mask_bdry(i,j-1) = 1.
!         endif
!       elseif (G%jdg_offset+j == gjec) then !top boundary
!         if (len_stress == 0. .OR. G%geoLonCv(i,j-1) <= len_stress) then
!           v_face_mask_bdry(i,j) = 0.
!         else
!           v_face_mask_bdry(i,j) = 1.
!         endif
!       endif

!       ! downstream boundary - CFBC

!       if (i+G%idg_offset == giec) then
!         u_face_mask_bdry(i,j) = 2.0
!       endif

!     enddo
!   enddo

!END MJH end subroutine initialize_ice_shelf_boundary_channel

end module MOM_ice_shelf_initialize
