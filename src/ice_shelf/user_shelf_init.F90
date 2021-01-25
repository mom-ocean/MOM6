!> This module specifies the initial values and evolving properties of the
!! MOM6 ice shelf, using user-provided code.
module user_shelf_init

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_time_manager,  only : time_type, set_time, time_type_to_real
use MOM_unit_scaling,  only : unit_scale_type

implicit none ; private

#include <MOM_memory.h>

public USER_initialize_shelf_mass, USER_update_shelf_mass
public USER_init_ice_thickness

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure for the user_ice_shelf module
type, public :: user_ice_shelf_CS ; private
  real :: Rho_ocean  !< The ocean's typical density [R ~> kg m-3].
  real :: max_draft  !< The maximum ocean draft of the ice shelf [Z ~> m].
  real :: min_draft  !< The minimum ocean draft of the ice shelf [Z ~> m].
  real :: flat_shelf_width !< The range over which the shelf is min_draft thick [km].
  real :: shelf_slope_scale !< The range over which the shelf slopes [km].
  real :: pos_shelf_edge_0 !< The x-position of the shelf edge at time 0 [km].
  real :: shelf_speed  !< The ice shelf speed of translation [km day-1]
  logical :: first_call = .true. !< If true, this module has not been called before.
end type user_ice_shelf_CS

contains

!> This subroutine sets up the initial mass and area covered by the ice shelf, based on user-provided code.
subroutine USER_initialize_shelf_mass(mass_shelf, area_shelf_h, h_shelf, hmask, G, US, CS, param_file, new_sim)

  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: mass_shelf !< The ice shelf mass per unit area averaged
                                                  !! over the full ocean cell [R Z ~> kg m-2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: hmask !< A mask indicating which tracer points are
                                                !! partly or fully covered by an ice-shelf
  type(unit_scale_type),   intent(in)  :: US    !< A structure containing unit conversion factors
  type(user_ice_shelf_CS), pointer     :: CS    !< A pointer to the user ice shelf control structure
  type(param_file_type),   intent(in)  :: param_file !< A structure to parse for run-time parameters
  logical,                 intent(in)  :: new_sim  !< If true, this is a new run; otherwise it is
                                                   !! being started from a restart file.

! This subroutine sets up the initial mass and area covered by the ice shelf.
  real :: max_draft  ! The maximum ocean draft of the ice shelf [Z ~> m].
  real :: min_draft  ! The minimum ocean draft of the ice shelf [Z ~> m].
  real :: flat_shelf_width ! The range over which the shelf is min_draft thick.
  real :: c1 ! The maximum depths in m.
  character(len=40) :: mdl = "USER_initialize_shelf_mass" ! This subroutine's name.
  integer :: i, j

  ! call MOM_error(FATAL, "USER_shelf_init.F90, USER_set_shelf_mass: " // &
  !  "Unmodified user routine called - you must edit the routine to use it")

  if (.not.associated(CS)) allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  if (CS%first_call) call write_user_log(param_file)
  CS%first_call = .false.
  call get_param(param_file, mdl, "RHO_0", CS%Rho_ocean, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINSEQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "SHELF_MAX_DRAFT", CS%max_draft, &
                 units="m", default=1.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "SHELF_MIN_DRAFT", CS%min_draft, &
                 units="m", default=1.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "FLAT_SHELF_WIDTH", CS%flat_shelf_width, &
                 units="axis_units", default=0.0)
  call get_param(param_file, mdl, "SHELF_SLOPE_SCALE", CS%shelf_slope_scale, &
                 units="axis_units", default=0.0)
  call get_param(param_file, mdl, "SHELF_EDGE_POS_0", CS%pos_shelf_edge_0, &
                 units="axis_units", default=0.0)
  call get_param(param_file, mdl, "SHELF_SPEED", CS%shelf_speed, &
                 units="axis_units day-1", default=0.0)

  call USER_update_shelf_mass(mass_shelf, area_shelf_h, h_shelf, hmask, G, CS, set_time(0,0), new_sim)

end subroutine USER_initialize_shelf_mass

!> This subroutine updates the ice shelf thickness, as specified by user-provided code.
subroutine USER_init_ice_thickness(h_shelf, area_shelf_h, hmask, G, US, param_file)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: h_shelf !< The ice shelf thickness [m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: hmask !< A mask indicating which tracer points are
                                                !! partly or fully covered by an ice-shelf
  type(unit_scale_type),   intent(in)  :: US    !< A structure containing unit conversion factors
  type(param_file_type),   intent(in)  :: param_file !< A structure to parse for run-time parameters

  ! This subroutine initializes the ice shelf thickness.  Currently it does so
  ! calling USER_initialize_shelf_mass, but this can be revised as needed.
  real, dimension(SZI_(G),SZJ_(G)) :: mass_shelf
  type(user_ice_shelf_CS), pointer :: CS => NULL()

  call USER_initialize_shelf_mass(mass_shelf, area_shelf_h, h_shelf, hmask, G, US, CS, param_file, .true.)

end subroutine USER_init_ice_thickness

!> This subroutine updates the ice shelf mass, as specified by user-provided code.
subroutine USER_update_shelf_mass(mass_shelf, area_shelf_h, h_shelf, hmask, G, CS, Time, new_sim)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(inout) :: mass_shelf !< The ice shelf mass per unit area averaged
                                                  !! over the full ocean cell [R Z ~> kg m-2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(inout) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(inout) :: hmask !< A mask indicating which tracer points are
                                                !! partly or fully covered by an ice-shelf
  type(user_ice_shelf_CS), pointer       :: CS   !< A pointer to the user ice shelf control structure
  type(time_type),         intent(in)    :: Time !< The current model time
  logical,                 intent(in)    :: new_sim !< If true, this the start of a new run.


  real :: c1, edge_pos, slope_pos
  integer :: i, j

  edge_pos = CS%pos_shelf_edge_0 + CS%shelf_speed*(time_type_to_real(Time) / 86400.0)

  slope_pos = edge_pos - CS%flat_shelf_width
  c1 = 0.0 ; if (CS%shelf_slope_scale > 0.0) c1 = 1.0 / CS%shelf_slope_scale


  do j=G%jsd,G%jed

    if (((j+G%jdg_offset) <= G%domain%njglobal+G%domain%njhalo) .AND. &
        ((j+G%jdg_offset) >= G%domain%njhalo+1)) then

      do i=G%isc,G%iec

!    if (((i+G%idg_offset) <= G%domain%niglobal+G%domain%nihalo) .AND. &
!           ((i+G%idg_offset) >= G%domain%nihalo+1)) then

        if ((j >= G%jsc) .and. (j <= G%jec)) then
          if (new_sim) then ; if (G%geoLonCu(i-1,j) >= edge_pos) then
            ! Everything past the edge is open ocean.
            mass_shelf(i,j) = 0.0
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
              h_shelf (i,j) = CS%min_draft
              mass_shelf(i,j) = CS%Rho_ocean * CS%min_draft
            else
              mass_shelf(i,j) = CS%Rho_ocean * (CS%min_draft + &
                     (CS%max_draft - CS%min_draft) * &
                     min(1.0, (c1*(slope_pos - G%geoLonT(i,j)))**2) )
              h_shelf(i,j) = (CS%min_draft + &
                     (CS%max_draft - CS%min_draft) * &
                     min(1.0, (c1*(slope_pos - G%geoLonT(i,j)))**2) )
            endif
          endif ; endif
        endif

        if ((i+G%idg_offset) == G%domain%nihalo+1) then
          hmask(i-1,j) = 3.0
        endif

      enddo
    endif
  enddo

end subroutine USER_update_shelf_mass

!> This subroutine writes out the user ice shelf code version number to the model log.
subroutine write_user_log(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters

  character(len=128) :: version = '$Id: user_shelf_init.F90,v 1.1.2.7 2012/06/19 22:15:52 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: MOM_ogrp $'
  character(len=40)  :: mdl = "user_shelf_init" ! This module's name.

  call log_version(param_file, mdl, version, tagname)

end subroutine write_user_log

end module user_shelf_init
