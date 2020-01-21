!> Reports integrated quantities for monitoring the model state
module MOM_sum_output

! This file is part of MOM6. See LICENSE.md for the license.

use iso_fortran_env, only : int64
use MOM_coms, only : sum_across_PEs, PE_here, root_PE, num_PEs, max_across_PEs
use MOM_coms, only : reproducing_sum, EFP_to_real, real_to_EFP
use MOM_coms, only : EFP_type, operator(+), operator(-), assignment(=)
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe, MOM_mesg
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_interface_heights, only : find_eta
use MOM_io, only : create_file, fieldtype, flush_file, open_file, reopen_file
use MOM_io, only : file_exists, slasher, vardesc, var_desc, write_field, get_filename_appendix
use MOM_io, only : APPEND_FILE, ASCII_FILE, SINGLE_FILE, WRITEONLY_FILE
use MOM_open_boundary, only : ocean_OBC_type, OBC_segment_type
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_time_manager, only : time_type, get_time, get_date, set_time, operator(>)
use MOM_time_manager, only : operator(+), operator(-), operator(*), operator(/)
use MOM_time_manager, only : operator(/=), operator(<=), operator(>=), operator(<)
use MOM_time_manager, only : get_calendar_type, time_type_to_real, NO_CALENDAR
use MOM_tracer_flow_control, only : tracer_flow_control_CS, call_tracer_stocks
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface, thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use mpp_mod, only : mpp_chksum

use netcdf

implicit none ; private

#include <MOM_memory.h>

public write_energy, accumulate_net_input, MOM_sum_output_init

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

integer, parameter :: NUM_FIELDS = 17 !< Number of diagnostic fields
character (*), parameter :: depth_chksum_attr = "bathyT_checksum"
                                      !< Checksum attribute name of G%bathyT
                                      !! over the compute domain
character (*), parameter :: area_chksum_attr = "mask2dT_areaT_checksum"
                                      !< Checksum attribute of name of
                                      !! G%mask2dT * G%areaT over the compute
                                      !! domain

!> A list of depths and corresponding globally integrated ocean area at each
!! depth and the ocean volume below each depth.
type :: Depth_List
  real :: depth       !< A depth [m].
  real :: area        !< The cross-sectional area of the ocean at that depth [m2].
  real :: vol_below   !< The ocean volume below that depth [m3].
end type Depth_List

!> The control structure for the MOM_sum_output module
type, public :: sum_output_CS ; private
  type(Depth_List), pointer, dimension(:) :: DL => NULL() !< The sorted depth list.
  integer :: list_size          !< length of sorting vector <= niglobal*njglobal

  integer, allocatable, dimension(:) :: lH
                                !< This saves the entry in DL with a volume just
                                !! less than the volume of fluid below the interface.
  logical :: do_APE_calc        !<   If true, calculate the available potential energy of the
                                !! interfaces.  Disabling this reduces the memory footprint of
                                !! high-PE-count models dramatically.
  logical :: read_depth_list    !<   Read the depth list from a file if it exists
                                !! and write it if it doesn't.
  character(len=200) :: depth_list_file  !< The name of the depth list file.
  real    :: D_list_min_inc     !<  The minimum increment [Z ~> m], between the depths of the
                                !! entries in the depth-list file, 0 by default.
  logical :: require_depth_list_chksum
                                !< Require matching checksums in Depth_list.nc when reading
                                !! the file.
  logical :: update_depth_list_chksum
                                !< Automatically update the Depth_list.nc file if the
                                !! checksums are missing or do not match current values.
  logical :: use_temperature    !<   If true, temperature and salinity are state variables.
  real    :: fresh_water_input  !<   The total mass of fresh water added by surface fluxes
                                !! since the last time that write_energy was called [kg].
  real    :: mass_prev          !<   The total ocean mass the last time that
                                !! write_energy was called [kg].
  real    :: salt_prev          !<   The total amount of salt in the ocean the last
                                !! time that write_energy was called [ppt kg].
  real    :: net_salt_input     !<   The total salt added by surface fluxes since the last
                                !! time that write_energy was called [ppt kg].
  real    :: heat_prev          !<  The total amount of heat in the ocean the last
                                !! time that write_energy was called [J].
  real    :: net_heat_input     !<  The total heat added by surface fluxes since the last
                                !! the last time that write_energy was called [J].
  type(EFP_type) :: fresh_water_in_EFP !< An extended fixed point version of fresh_water_input
  type(EFP_type) :: net_salt_in_EFP !< An extended fixed point version of net_salt_input
  type(EFP_type) :: net_heat_in_EFP !< An extended fixed point version of net_heat_input
  type(EFP_type) :: heat_prev_EFP !< An extended fixed point version of heat_prev
  type(EFP_type) :: salt_prev_EFP !< An extended fixed point version of salt_prev
  type(EFP_type) :: mass_prev_EFP !< An extended fixed point version of mass_prev
  real    :: dt_in_T            !< The baroclinic dynamics time step [T ~> s].

  type(time_type) :: energysavedays            !< The interval between writing the energies
                                               !! and other integral quantities of the run.
  type(time_type) :: energysavedays_geometric  !< The starting interval for computing a geometric
                                               !! progression of time deltas between calls to
                                               !! write_energy. This interval will increase by a factor of 2.
                                               !! after each call to write_energy.
  logical         :: energysave_geometric      !< Logical to control whether calls to write_energy should
                                               !! follow a geometric progression
  type(time_type) :: write_energy_time         !< The next time to write to the energy file.
  type(time_type) :: geometric_end_time        !< Time at which to stop the geometric progression
                                               !! of calls to write_energy and revert to the standard
                                               !! energysavedays interval

  real    :: timeunit           !<  The length of the units for the time axis [s].
  logical :: date_stamped_output !< If true, use dates (not times) in messages to stdout.
  type(time_type) :: Start_time !< The start time of the simulation.
                                ! Start_time is set in MOM_initialization.F90
  integer, pointer :: ntrunc => NULL() !< The number of times the velocity has been
                                !! truncated since the last call to write_energy.
  real    :: max_Energy         !< The maximum permitted energy per unit mass.  If there is
                                !! more energy than this, the model should stop [m2 s-2].
  integer :: maxtrunc           !< The number of truncations per energy save
                                !! interval at which the run is stopped.
  logical :: write_stocks       !< If true, write the integrated tracer amounts
                                !! to stdout when the energy files are written.
  integer :: previous_calls = 0 !< The number of times write_energy has been called.
  integer :: prev_n = 0         !< The value of n from the last call.
  integer :: fileenergy_nc      !< NetCDF id of the energy file.
  integer :: fileenergy_ascii   !< The unit number of the ascii version of the energy file.
  type(fieldtype), dimension(NUM_FIELDS+MAX_FIELDS_) :: &
             fields             !< fieldtype variables for the output fields.
  character(len=200) :: energyfile  !< The name of the energy file with path.
end type sum_output_CS

contains

!> MOM_sum_output_init initializes the parameters and settings for the MOM_sum_output module.
subroutine MOM_sum_output_init(G, US, param_file, directory, ntrnc, &
                               Input_start_time, CS)
  type(ocean_grid_type),  intent(in)    :: G          !< The ocean's grid structure.
  type(unit_scale_type),  intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time
                                                      !! parameters.
  character(len=*),       intent(in)    :: directory  !< The directory where the energy file goes.
  integer, target,        intent(inout) :: ntrnc      !< The integer that stores the number of times
                                                      !! the velocity has been truncated since the
                                                      !! last call to write_energy.
  type(time_type),        intent(in)    :: Input_start_time !< The start time of the simulation.
  type(Sum_output_CS),    pointer       :: CS         !< A pointer that is set to point to the
                                                      !! control structure for this module.
  ! Local variables
  real :: Time_unit ! The time unit in seconds for ENERGYSAVEDAYS.
  real :: Rho_0     ! A reference density [kg m-3]
  real :: maxvel    ! The maximum permitted velocity [m s-1]
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_sum_output" ! This module's name.
  character(len=200) :: energyfile  ! The name of the energy file.
  character(len=32) :: filename_appendix = '' !fms appendix to filename for ensemble runs

  if (associated(CS)) then
    call MOM_error(WARNING, "MOM_sum_output_init called with associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "CALCULATE_APE", CS%do_APE_calc, &
                 "If true, calculate the available potential energy of "//&
                 "the interfaces.  Setting this to false reduces the "//&
                 "memory footprint of high-PE-count models dramatically.", &
                 default=.true.)
  call get_param(param_file, mdl, "WRITE_STOCKS", CS%write_stocks, &
                 "If true, write the integrated tracer amounts to stdout "//&
                 "when the energy files are written.", default=.true.)
  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state "//&
                 "variables.", default=.true.)
  call get_param(param_file, mdl, "DT", CS%dt_in_T, &
                 "The (baroclinic) dynamics time step.", &
                 units="s", scale=US%s_to_T, fail_if_missing=.true.)
  call get_param(param_file, mdl, "MAXTRUNC", CS%maxtrunc, &
                 "The run will be stopped, and the day set to a very "//&
                 "large value if the velocity is truncated more than "//&
                 "MAXTRUNC times between energy saves.  Set MAXTRUNC to 0 "//&
                 "to stop if there is any truncation of velocities.", &
                 units="truncations save_interval-1", default=0)

  call get_param(param_file, mdl, "MAX_ENERGY", CS%max_Energy, &
                 "The maximum permitted average energy per unit mass; the "//&
                 "model will be stopped if there is more energy than "//&
                 "this.  If zero or negative, this is set to 10*MAXVEL^2.", &
                 units="m2 s-2", default=0.0)
  if (CS%max_Energy <= 0.0) then
    call get_param(param_file, mdl, "MAXVEL", maxvel, &
                 "The maximum velocity allowed before the velocity "//&
                 "components are truncated.", units="m s-1", default=3.0e8)
    CS%max_Energy = 10.0 * maxvel**2
    call log_param(param_file, mdl, "MAX_ENERGY as used", CS%max_Energy)
  endif

  call get_param(param_file, mdl, "ENERGYFILE", energyfile, &
                 "The file to use to write the energies and globally "//&
                 "summed diagnostics.", default="ocean.stats")

  !query fms_io if there is a filename_appendix (for ensemble runs)
  call get_filename_appendix(filename_appendix)
  if (len_trim(filename_appendix) > 0) then
     energyfile = trim(energyfile) //'.'//trim(filename_appendix)
  endif

  CS%energyfile = trim(slasher(directory))//trim(energyfile)
  call log_param(param_file, mdl, "output_path/ENERGYFILE", CS%energyfile)
#ifdef STATSLABEL
  CS%energyfile = trim(CS%energyfile)//"."//trim(adjustl(STATSLABEL))
#endif

  call get_param(param_file, mdl, "DATE_STAMPED_STDOUT", CS%date_stamped_output, &
                 "If true, use dates (not times) in messages to stdout", &
                 default=.true.)
  call get_param(param_file, mdl, "TIMEUNIT", CS%Timeunit, &
                 "The time unit in seconds a number of input fields", &
                 units="s", default=86400.0)
  if (CS%Timeunit < 0.0) CS%Timeunit = 86400.0



  if (CS%do_APE_calc) then
    call get_param(param_file, mdl, "READ_DEPTH_LIST", CS%read_depth_list, &
                   "Read the depth list from a file if it exists or "//&
                   "create that file otherwise.", default=.false.)
    call get_param(param_file, mdl, "DEPTH_LIST_MIN_INC", CS%D_list_min_inc, &
                   "The minimum increment between the depths of the "//&
                   "entries in the depth-list file.", &
                   units="m", default=1.0E-10, scale=US%m_to_Z)
    if (CS%read_depth_list) then
      call get_param(param_file, mdl, "DEPTH_LIST_FILE", CS%depth_list_file, &
                   "The name of the depth list file.", default="Depth_list.nc")
      if (scan(CS%depth_list_file,'/') == 0) &
        CS%depth_list_file = trim(slasher(directory))//trim(CS%depth_list_file)

      call get_param(param_file, mdl, "REQUIRE_DEPTH_LIST_CHECKSUMS", &
                     CS%require_depth_list_chksum, &
                 "Require that matching checksums be in Depth_list.nc "//&
                 "when reading the file.", default=.true.)
      if (.not. CS%require_depth_list_chksum) &
        call get_param(param_file, mdl, "UPDATE_DEPTH_LIST_CHECKSUMS", &
                     CS%update_depth_list_chksum, &
                 "Automatically update the Depth_list.nc file if the "//&
                 "checksums are missing or do not match current values.", &
                 default=.false.)
    endif

    allocate(CS%lH(G%ke))
    call depth_list_setup(G, US, CS)
  else
    CS%list_size = 0
  endif

  call get_param(param_file, mdl, "TIMEUNIT", Time_unit, &
                 "The time unit for ENERGYSAVEDAYS.", &
                 units="s", default=86400.0)
  call get_param(param_file, mdl, "ENERGYSAVEDAYS",CS%energysavedays, &
                 "The interval in units of TIMEUNIT between saves of the "//&
                 "energies of the run and other globally summed diagnostics.",&
                 default=set_time(0,days=1), timeunit=Time_unit)
  call get_param(param_file, mdl, "ENERGYSAVEDAYS_GEOMETRIC",CS%energysavedays_geometric, &
                 "The starting interval in units of TIMEUNIT for the first call "//&
                 "to save the energies of the run and other globally summed diagnostics. "//&
                 "The interval increases by a factor of 2. after each call to write_energy.",&
                 default=set_time(seconds=0), timeunit=Time_unit)

  if ((time_type_to_real(CS%energysavedays_geometric) > 0.) .and. &
     (CS%energysavedays_geometric < CS%energysavedays)) then
         CS%energysave_geometric = .true.
  else
         CS%energysave_geometric = .false.
  endif

  CS%Start_time = Input_start_time
  CS%ntrunc => ntrnc

end subroutine MOM_sum_output_init

!> MOM_sum_output_end deallocates memory used by the MOM_sum_output module.
subroutine MOM_sum_output_end(CS)
  type(Sum_output_CS), pointer :: CS  !< The control structure returned by a
                                      !! previous call to MOM_sum_output_init.
  if (associated(CS)) then
    if (CS%do_APE_calc) then
      deallocate(CS%lH, CS%DL)
    endif

    deallocate(CS)
  endif
end subroutine MOM_sum_output_end

!>  This subroutine calculates and writes the total model energy, the energy and
!! mass of each layer, and other globally integrated  physical quantities.
subroutine write_energy(u, v, h, tv, day, n, G, GV, US, CS, tracer_CSp, OBC, dt_forcing)
  type(ocean_grid_type),   intent(in)    :: G   !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: u   !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(in)    :: v   !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                           intent(in)    :: h   !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv  !< A structure pointing to various
                                                !! thermodynamic variables.
  type(time_type),         intent(in)    :: day !< The current model time.
  integer,                 intent(in)    :: n   !< The time step number of the
                                                !! current execution.
  type(Sum_output_CS),     pointer       :: CS  !< The control structure returned by a
                                                !! previous call to MOM_sum_output_init.
  type(tracer_flow_control_CS), &
                    optional, pointer    :: tracer_CSp !< tracer control structure.
  type(ocean_OBC_type),         &
                    optional, pointer    :: OBC !< Open boundaries control structure.
  type(time_type),  optional, intent(in) :: dt_forcing !< The forcing time step
  ! Local variables
  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! The height of interfaces [Z ~> m].
  real :: areaTm(SZI_(G),SZJ_(G)) ! A masked version of areaT [m2].
  real :: KE(SZK_(G))  ! The total kinetic energy of a layer [J].
  real :: PE(SZK_(G)+1)! The available potential energy of an interface [J].
  real :: KE_tot       ! The total kinetic energy [J].
  real :: PE_tot       ! The total available potential energy [J].
  real :: Z_0APE(SZK_(G)+1) ! The uniform depth which overlies the same
                       ! volume as is below an interface [Z ~> m].
  real :: H_0APE(SZK_(G)+1) ! A version of Z_0APE, converted to m, usually positive.
  real :: toten        ! The total kinetic & potential energies of
                       ! all layers [J] (i.e. kg m2 s-2).
  real :: En_mass      ! The total kinetic and potential energies divided by
                       ! the total mass of the ocean [m2 s-2].
  real :: vol_lay(SZK_(G))  ! The volume of fluid in a layer [Z m2 ~> m3].
  real :: volbelow     ! The volume of all layers beneath an interface [Z m2 ~> m3].
  real :: mass_lay(SZK_(G)) ! The mass of fluid in a layer [kg].
  real :: mass_tot     ! The total mass of the ocean [kg].
  real :: vol_tot      ! The total ocean volume [m3].
  real :: mass_chg     ! The change in total ocean mass of fresh water since
                       ! the last call to this subroutine [kg].
  real :: mass_anom    ! The change in fresh water that cannot be accounted for
                       ! by the surface fluxes [kg].
  real :: Salt         ! The total amount of salt in the ocean [ppt kg].
  real :: Salt_chg     ! The change in total ocean salt since the last call
                       ! to this subroutine [ppt kg].
  real :: Salt_anom    ! The change in salt that cannot be accounted for by
                       ! the surface fluxes [ppt kg].
  real :: salin        ! The mean salinity of the ocean [ppt].
  real :: salin_chg    ! The change in total salt since the last call
                       ! to this subroutine divided by total mass [ppt].
  real :: salin_anom   ! The change in total salt that cannot be accounted for by
                       ! the surface fluxes divided by total mass [ppt].
  real :: salin_mass_in ! The mass of salt input since the last call [kg].
  real :: Heat         ! The total amount of Heat in the ocean [J].
  real :: Heat_chg     ! The change in total ocean heat since the last call
                       ! to this subroutine [J].
  real :: Heat_anom    ! The change in heat that cannot be accounted for by
                       ! the surface fluxes [J].
  real :: temp         ! The mean potential temperature of the ocean [degC].
  real :: temp_chg     ! The change in total heat divided by total heat capacity
                       ! of the ocean since the last call to this subroutine, degC.
  real :: temp_anom    ! The change in total heat that cannot be accounted for
                       ! by the surface fluxes, divided by the total heat
                       ! capacity of the ocean [degC].
  real :: hint         ! The deviation of an interface from H [Z ~> m].
  real :: hbot         ! 0 if the basin is deeper than H, or the
                       ! height of the basin depth over H otherwise [Z ~> m].
                       ! This makes PE only include real fluid.
  real :: hbelow       ! The depth of fluid in all layers beneath an interface [Z ~> m].
  type(EFP_type) :: &
    mass_EFP, &        ! Extended fixed point sums of total mass, etc.
    salt_EFP, heat_EFP, salt_chg_EFP, heat_chg_EFP, mass_chg_EFP, &
    mass_anom_EFP, salt_anom_EFP, heat_anom_EFP
  real :: CFL_trans    ! A transport-based definition of the CFL number [nondim].
  real :: CFL_lin      ! A simpler definition of the CFL number [nondim].
  real :: max_CFL(2)   ! The maxima of the CFL numbers [nondim].
  real :: Irho0        ! The inverse of the reference density [m3 kg-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    tmp1               ! A temporary array
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: &
    PE_pt              ! The potential energy at each point [J].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    Temp_int, Salt_int ! Layer and cell integrated heat and salt [J] and [g Salt].
  real :: H_to_kg_m2   ! Local copy of a unit conversion factor.
  real :: KE_scale_factor   ! The combination of unit rescaling factors in the kinetic energy
                            ! calculation [kg T2 L-2 s-2 H-1 ~> kg m-3 or nondim]
  integer :: num_nc_fields  ! The number of fields that will actually go into
                            ! the NetCDF file.
  integer :: i, j, k, is, ie, js, je, ns, nz, m, Isq, Ieq, Jsq, Jeq
  integer :: l, lbelow, labove   ! indices of deep_area_vol, used to find Z_0APE.
                                 ! lbelow & labove are lower & upper limits for l
                                 ! in the search for the entry in lH to use.
  integer :: start_of_day, num_days
  real    :: reday, var
  character(len=240) :: energypath_nc
  character(len=200) :: mesg
  character(len=32)  :: mesg_intro, time_units, day_str, n_str, date_str
  logical :: date_stamped
  type(time_type) :: dt_force ! A time_type version of the forcing timestep.
  real :: Tr_stocks(MAX_FIELDS_)
  real :: Tr_min(MAX_FIELDS_), Tr_max(MAX_FIELDS_)
  real :: Tr_min_x(MAX_FIELDS_), Tr_min_y(MAX_FIELDS_), Tr_min_z(MAX_FIELDS_)
  real :: Tr_max_x(MAX_FIELDS_), Tr_max_y(MAX_FIELDS_), Tr_max_z(MAX_FIELDS_)
  logical :: Tr_minmax_got(MAX_FIELDS_) = .false.
  character(len=40), dimension(MAX_FIELDS_) :: &
    Tr_names, Tr_units
  integer :: nTr_stocks
  integer :: iyear, imonth, iday, ihour, iminute, isecond, itick ! For call to get_date()
  logical :: local_open_BC
  type(OBC_segment_type), pointer :: segment => NULL()

 ! A description for output of each of the fields.
  type(vardesc) :: vars(NUM_FIELDS+MAX_FIELDS_)

  ! write_energy_time is the next integral multiple of energysavedays.
  dt_force = set_time(seconds=2) ; if (present(dt_forcing)) dt_force = dt_forcing
  if (CS%previous_calls == 0) then
    if (CS%energysave_geometric) then
      if (CS%energysavedays_geometric < CS%energysavedays) then
        CS%write_energy_time = day + CS%energysavedays_geometric
        CS%geometric_end_time = CS%Start_time + CS%energysavedays * &
          (1 + (day - CS%Start_time) / CS%energysavedays)
      else
        CS%write_energy_time = CS%Start_time + CS%energysavedays * &
          (1 + (day - CS%Start_time) / CS%energysavedays)
      endif
    else
      CS%write_energy_time = CS%Start_time + CS%energysavedays * &
        (1 + (day - CS%Start_time) / CS%energysavedays)
    endif
  elseif (day + (dt_force/2) <= CS%write_energy_time) then
    return  ! Do not write this step
  else ! Determine the next write time before proceeding
    if (CS%energysave_geometric) then
      if (CS%write_energy_time + CS%energysavedays_geometric >= &
          CS%geometric_end_time) then
        CS%write_energy_time = CS%geometric_end_time
        CS%energysave_geometric = .false.  ! stop geometric progression
      else
        CS%write_energy_time = CS%write_energy_time + CS%energysavedays_geometric
      endif
      CS%energysavedays_geometric = CS%energysavedays_geometric*2
    else
      CS%write_energy_time = CS%write_energy_time + CS%energysavedays
    endif
  endif

  num_nc_fields = 17
  if (.not.CS%use_temperature) num_nc_fields = 11
  vars(1) = var_desc("Ntrunc","Nondim","Number of Velocity Truncations",'1','1')
  vars(2) = var_desc("En","Joules","Total Energy",'1','1')
  vars(3) = var_desc("APE","Joules","Total Interface APE",'1','i')
  vars(4) = var_desc("KE","Joules","Total Layer KE",'1','L')
  vars(5) = var_desc("H0","meter","Zero APE Depth of Interface",'1','i')
  vars(6) = var_desc("Mass_lay","kg","Total Layer Mass",'1','L')
  vars(7) = var_desc("Mass","kg","Total Mass",'1','1')
  vars(8) = var_desc("Mass_chg","kg","Total Mass Change between Entries",'1','1')
  vars(9) = var_desc("Mass_anom","kg","Anomalous Total Mass Change",'1','1')
  vars(10) = var_desc("max_CFL_trans","Nondim","Maximum finite-volume CFL",'1','1')
  vars(11) = var_desc("max_CFL_lin","Nondim","Maximum finite-difference CFL",'1','1')
  if (CS%use_temperature) then
    vars(12) = var_desc("Salt","kg","Total Salt",'1','1')
    vars(13) = var_desc("Salt_chg","kg","Total Salt Change between Entries",'1','1')
    vars(14) = var_desc("Salt_anom","kg","Anomalous Total Salt Change",'1','1')
    vars(15) = var_desc("Heat","Joules","Total Heat",'1','1')
    vars(16) = var_desc("Heat_chg","Joules","Total Heat Change between Entries",'1','1')
    vars(17) = var_desc("Heat_anom","Joules","Anomalous Total Heat Change",'1','1')
  endif

  local_open_BC = .false.
  if (present(OBC)) then ; if (associated(OBC)) then
    local_open_BC = (OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally)
  endif ; endif

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  H_to_kg_m2 = GV%H_to_kg_m2

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "write_energy: Module must be initialized before it is used.")

  do j=js,je ; do i=is,ie
    areaTm(i,j) = G%mask2dT(i,j)*US%L_to_m**2*G%areaT(i,j)
  enddo ; enddo

  if (GV%Boussinesq) then
    tmp1(:,:,:) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      tmp1(i,j,k) = h(i,j,k) * (H_to_kg_m2*areaTm(i,j))
    enddo ; enddo ; enddo

    ! This block avoids using the points beyond an open boundary condition
    ! in the accumulation of mass, but perhaps it would be unnecessary if there
    ! were a more judicious use of masks in the loops 4 or 7 lines above.
    if (local_open_BC) then
      do ns=1, OBC%number_of_segments
        segment => OBC%segment(ns)
        if (.not. segment%on_pe .or. segment%specified) cycle
        I=segment%HI%IsdB ; J=segment%HI%JsdB
        if (segment%direction == OBC_DIRECTION_E) then
          do k=1,nz ; do j=segment%HI%jsd,segment%HI%jed
            tmp1(i+1,j,k) = 0.0
          enddo ; enddo
        elseif (segment%direction == OBC_DIRECTION_W) then
          do k=1,nz ; do j=segment%HI%jsd,segment%HI%jed
            tmp1(i,j,k) = 0.0
          enddo ; enddo
        elseif (segment%direction == OBC_DIRECTION_N) then
          do k=1,nz ; do i=segment%HI%isd,segment%HI%ied
            tmp1(i,j+1,k) = 0.0
          enddo ; enddo
        elseif (segment%direction == OBC_DIRECTION_S) then
          do k=1,nz ; do i=segment%HI%isd,segment%HI%ied
            tmp1(i,j,k) = 0.0
          enddo ; enddo
        endif
      enddo
    endif

    mass_tot = reproducing_sum(tmp1, sums=mass_lay, EFP_sum=mass_EFP)
    do k=1,nz ; vol_lay(k) = (GV%H_to_Z/H_to_kg_m2)*mass_lay(k) ; enddo
  else
    tmp1(:,:,:) = 0.0
    if (CS%do_APE_calc) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        tmp1(i,j,k) = H_to_kg_m2 * h(i,j,k) * areaTm(i,j)
      enddo ; enddo ; enddo
      mass_tot = reproducing_sum(tmp1, sums=mass_lay, EFP_sum=mass_EFP)

      call find_eta(h, tv, G, GV, US, eta)
      do k=1,nz ; do j=js,je ; do i=is,ie
        tmp1(i,j,k) = US%Z_to_m*(eta(i,j,K)-eta(i,j,K+1)) * areaTm(i,j)
      enddo ; enddo ; enddo
      vol_tot = reproducing_sum(tmp1, sums=vol_lay)
      do k=1,nz ; vol_lay(k) = US%m_to_Z * vol_lay(k) ; enddo
    else
      do k=1,nz ; do j=js,je ; do i=is,ie
        tmp1(i,j,k) = H_to_kg_m2 * h(i,j,k) * areaTm(i,j)
      enddo ; enddo ; enddo
      mass_tot = reproducing_sum(tmp1, sums=mass_lay, EFP_sum=mass_EFP)
      do k=1,nz ; vol_lay(k) = US%m_to_Z * (mass_lay(k) / (US%R_to_kg_m3*GV%Rho0)) ; enddo
    endif
  endif ! Boussinesq

  nTr_stocks = 0
  if (present(tracer_CSp)) then
    call call_tracer_stocks(h, Tr_stocks, G, GV, tracer_CSp, stock_names=Tr_names, &
                            stock_units=Tr_units, num_stocks=nTr_stocks,&
                            got_min_max=Tr_minmax_got, global_min=Tr_min, global_max=Tr_max, &
                            xgmin=Tr_min_x, ygmin=Tr_min_y, zgmin=Tr_min_z,&
                            xgmax=Tr_max_x, ygmax=Tr_max_y, zgmax=Tr_max_z)
    if (nTr_stocks > 0) then
      do m=1,nTr_stocks
        vars(num_nc_fields+m) = var_desc(Tr_names(m), units=Tr_units(m), &
                      longname=Tr_names(m), hor_grid='1', z_grid='1')
      enddo
      num_nc_fields = num_nc_fields + nTr_stocks
    endif
  endif

  if (CS%previous_calls == 0) then
    CS%mass_prev = mass_tot ; CS%fresh_water_input = 0.0

    CS%mass_prev_EFP = mass_EFP
    CS%fresh_water_in_EFP = real_to_EFP(0.0)

    !  Reopen or create a text output file, with an explanatory header line.
    if (is_root_pe()) then
      if (day > CS%Start_time) then
        call open_file(CS%fileenergy_ascii, trim(CS%energyfile), &
                       action=APPEND_FILE, form=ASCII_FILE, nohdrs=.true.)
      else
        call open_file(CS%fileenergy_ascii, trim(CS%energyfile), &
                       action=WRITEONLY_FILE, form=ASCII_FILE, nohdrs=.true.)
        if (abs(CS%timeunit - 86400.0) < 1.0) then
          if (CS%use_temperature) then
            write(CS%fileenergy_ascii,'("  Step,",7x,"Day,  Truncs,      &
                &Energy/Mass,      Maximum CFL,  Mean Sea Level,  Total Mass,  Mean Salin, &
                &Mean Temp, Frac Mass Err,   Salin Err,    Temp Err")')
            write(CS%fileenergy_ascii,'(12x,"[days]",17x,"[m2 s-2]",11x,"[Nondim]",7x,"[m]",13x,&
                &"[kg]",9x,"[PSU]",6x,"[degC]",7x,"[Nondim]",8x,"[PSU]",8x,"[degC]")')
          else
            write(CS%fileenergy_ascii,'("  Step,",7x,"Day,  Truncs,      &
                &Energy/Mass,      Maximum CFL,  Mean sea level,   Total Mass,    Frac Mass Err")')
            write(CS%fileenergy_ascii,'(12x,"[days]",17x,"[m2 s-2]",11x,"[Nondim]",8x,"[m]",13x,&
                &"[kg]",11x,"[Nondim]")')
          endif
        else
          if ((CS%timeunit >= 0.99) .and. (CS%timeunit < 1.01)) then
            time_units = "           [seconds]     "
          elseif ((CS%timeunit >= 3599.0) .and. (CS%timeunit < 3601.0)) then
            time_units = "            [hours]      "
          elseif ((CS%timeunit >= 86399.0) .and. (CS%timeunit < 86401.0)) then
            time_units = "             [days]      "
          elseif ((CS%timeunit >= 3.0e7) .and. (CS%timeunit < 3.2e7)) then
            time_units = "            [years]      "
          else
            write(time_units,'(9x,"[",es8.2," s]    ")') CS%timeunit
          endif

          if (CS%use_temperature) then
            write(CS%fileenergy_ascii,'("  Step,",7x,"Time, Truncs,      &
                &Energy/Mass,      Maximum CFL,  Mean Sea Level,  Total Mass,  Mean Salin, &
                &Mean Temp, Frac Mass Err,   Salin Err,    Temp Err")')
            write(CS%fileenergy_ascii,'(A25,10x,"[m2 s-2]",11x,"[Nondim]",7x,"[m]",13x,&
                &"[kg]",9x,"[PSU]",6x,"[degC]",7x,"[Nondim]",8x,"[PSU]",6x,&
                &"[degC]")') time_units
          else
            write(CS%fileenergy_ascii,'("  Step,",7x,"Time, Truncs,      &
                &Energy/Mass,      Maximum CFL,  Mean sea level,   Total Mass,    Frac Mass Err")')
            write(CS%fileenergy_ascii,'(A25,10x,"[m2 s-2]",11x,"[Nondim]",8x,"[m]",13x,&
                &"[kg]",11x,"[Nondim]")') time_units
          endif
        endif
      endif
    endif

    energypath_nc = trim(CS%energyfile) // ".nc"
    if (day > CS%Start_time) then
      call reopen_file(CS%fileenergy_nc, trim(energypath_nc), vars, &
                       num_nc_fields, CS%fields, SINGLE_FILE, CS%timeunit, &
                       G=G, GV=GV)
    else
      call create_file(CS%fileenergy_nc, trim(energypath_nc), vars, &
                       num_nc_fields, CS%fields, SINGLE_FILE, CS%timeunit, &
                       G=G, GV=GV)
    endif
  endif

  if (CS%do_APE_calc) then
    lbelow = 1 ; volbelow = 0.0
    do k=nz,1,-1
      volbelow = volbelow + vol_lay(k)
      if ((volbelow >= CS%DL(CS%lH(k))%vol_below) .and. &
          (volbelow < CS%DL(CS%lH(k)+1)%vol_below)) then
        l = CS%lH(k)
      else
        labove=CS%list_size+1
        l = (labove + lbelow) / 2
        do while (l > lbelow)
          if (volbelow < CS%DL(l)%vol_below) then ; labove = l
          else ; lbelow = l ; endif
          l = (labove + lbelow) / 2
        enddo
        CS%lH(k) = l
      endif
      lbelow = l
      Z_0APE(K) = CS%DL(l)%depth - (volbelow - CS%DL(l)%vol_below) / CS%DL(l)%area
    enddo
    Z_0APE(nz+1) = CS%DL(2)%depth

    !   Calculate the Available Potential Energy integrated over each
    ! interface.  With a nonlinear equation of state or with a bulk
    ! mixed layer this calculation is only approximate.  With an ALE model
    ! this does not make sense.
    PE_pt(:,:,:) = 0.0
    if (GV%Boussinesq) then
      do j=js,je ; do i=is,ie
        hbelow = 0.0
        do k=nz,1,-1
          hbelow = hbelow + h(i,j,k) * GV%H_to_Z
          hint = Z_0APE(K) + (hbelow - G%bathyT(i,j))
          hbot = Z_0APE(K) - G%bathyT(i,j)
          hbot = (hbot + ABS(hbot)) * 0.5
          PE_pt(i,j,K) = 0.5 * areaTm(i,j) * US%Z_to_m*US%L_T_to_m_s**2*(US%R_to_kg_m3*GV%Rho0*GV%g_prime(K)) * &
                  (hint * hint - hbot * hbot)
        enddo
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        do k=nz,1,-1
          hint = Z_0APE(K) + eta(i,j,K)  ! eta and H_0 have opposite signs.
          hbot = max(Z_0APE(K) - G%bathyT(i,j), 0.0)
          PE_pt(i,j,K) = 0.5 * (areaTm(i,j) * US%Z_to_m*US%L_T_to_m_s**2*(US%R_to_kg_m3*GV%Rho0*GV%g_prime(K))) * &
                  (hint * hint - hbot * hbot)
        enddo
      enddo ; enddo
    endif

    PE_tot = reproducing_sum(PE_pt, sums=PE)
    do k=1,nz+1 ; H_0APE(K) = US%Z_to_m*Z_0APE(K) ; enddo
  else
    PE_tot = 0.0
    do k=1,nz+1 ; PE(K) = 0.0 ; H_0APE(K) = 0.0 ; enddo
  endif

! Calculate the Kinetic Energy integrated over each layer.
  KE_scale_factor = GV%H_to_kg_m2*US%L_T_to_m_s**2
  tmp1(:,:,:) = 0.0
  do k=1,nz ; do j=js,je ; do i=is,ie
    tmp1(i,j,k) = (0.25 * KE_scale_factor * (areaTm(i,j) * h(i,j,k))) * &
            (u(I-1,j,k)**2 + u(I,j,k)**2 + v(i,J-1,k)**2 + v(i,J,k)**2)
  enddo ; enddo ; enddo
  KE_tot = reproducing_sum(tmp1, sums=KE)

  toten = KE_tot + PE_tot

  Salt = 0.0 ; Heat = 0.0
  if (CS%use_temperature) then
    Temp_int(:,:) = 0.0 ; Salt_int(:,:) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      Salt_int(i,j) = Salt_int(i,j) + tv%S(i,j,k) * &
                      (h(i,j,k)*(H_to_kg_m2 * areaTm(i,j)))
      Temp_int(i,j) = Temp_int(i,j) + (tv%C_p * tv%T(i,j,k)) * &
                      (h(i,j,k)*(H_to_kg_m2 * areaTm(i,j)))
    enddo ; enddo ; enddo
    Salt = reproducing_sum(Salt_int, EFP_sum=salt_EFP)
    Heat = reproducing_sum(Temp_int, EFP_sum=heat_EFP)
  endif

! Calculate the maximum CFL numbers.
  max_CFL(1:2) = 0.0
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    if (u(I,j,k) < 0.0) then
      CFL_trans = (-u(I,j,k) * CS%dt_in_T) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
    else
      CFL_trans = (u(I,j,k) * CS%dt_in_T) * (G%dy_Cu(I,j) * G%IareaT(i,j))
    endif
    CFL_lin = abs(u(I,j,k) * CS%dt_in_T) * G%IdxCu(I,j)
    max_CFL(1) = max(max_CFL(1), CFL_trans)
    max_CFL(2) = max(max_CFL(2), CFL_lin)
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    if (v(i,J,k) < 0.0) then
      CFL_trans = (-v(i,J,k) * CS%dt_in_T) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
    else
      CFL_trans = (v(i,J,k) * CS%dt_in_T) * (G%dx_Cv(i,J) * G%IareaT(i,j))
    endif
    CFL_lin = abs(v(i,J,k) * CS%dt_in_T) * G%IdyCv(i,J)
    max_CFL(1) = max(max_CFL(1), CFL_trans)
    max_CFL(2) = max(max_CFL(2), CFL_lin)
  enddo ; enddo ; enddo

  call sum_across_PEs(CS%ntrunc)
  !   Sum the various quantities across all the processors.  This sum is NOT
  ! guaranteed to be bitwise reproducible, even on the same decomposition.
  !   The sum of Tr_stocks should be reimplemented using the reproducing sums.
  if (nTr_stocks > 0) call sum_across_PEs(Tr_stocks,nTr_stocks)

  call max_across_PEs(max_CFL(1))
  call max_across_PEs(max_CFL(2))
  if (CS%use_temperature .and. CS%previous_calls == 0) then
    CS%salt_prev = Salt ; CS%net_salt_input = 0.0
    CS%heat_prev = Heat ; CS%net_heat_input = 0.0

    CS%salt_prev_EFP = salt_EFP ; CS%net_salt_in_EFP = real_to_EFP(0.0)
    CS%heat_prev_EFP = heat_EFP ; CS%net_heat_in_EFP = real_to_EFP(0.0)
  endif
  Irho0 = 1.0 / (US%R_to_kg_m3*GV%Rho0)

  if (CS%use_temperature) then
    Salt_chg_EFP = Salt_EFP - CS%salt_prev_EFP
    Salt_anom_EFP = Salt_chg_EFP - CS%net_salt_in_EFP
    Salt_chg = EFP_to_real(Salt_chg_EFP) ; Salt_anom = EFP_to_real(Salt_anom_EFP)
    Heat_chg_EFP = Heat_EFP - CS%heat_prev_EFP
    Heat_anom_EFP = Heat_chg_EFP - CS%net_heat_in_EFP
    Heat_chg = EFP_to_real(Heat_chg_EFP) ; Heat_anom = EFP_to_real(Heat_anom_EFP)
  endif

  mass_chg_EFP = mass_EFP - CS%mass_prev_EFP
  salin_mass_in = 0.0
  if (GV%Boussinesq) then
    mass_anom_EFP = mass_chg_EFP - CS%fresh_water_in_EFP
  else
    ! net_salt_input needs to be converted from ppt m s-1 to kg m-2 s-1.
    mass_anom_EFP = mass_chg_EFP - CS%fresh_water_in_EFP
    if (CS%use_temperature) &
      salin_mass_in = 0.001*EFP_to_real(CS%net_salt_in_EFP)
  endif
  mass_chg = EFP_to_real(mass_chg_EFP)
  mass_anom = EFP_to_real(mass_anom_EFP) - salin_mass_in

  if (CS%use_temperature) then
    salin = Salt / mass_tot ; salin_anom = Salt_anom / mass_tot
   ! salin_chg = Salt_chg / mass_tot
    temp = heat / (mass_tot*tv%C_p) ; temp_anom = Heat_anom / (mass_tot*tv%C_p)
  endif
  En_mass = toten / mass_tot

  call get_time(day, start_of_day, num_days)
  date_stamped = (CS%date_stamped_output .and. (get_calendar_type() /= NO_CALENDAR))
  if (date_stamped) &
    call get_date(day, iyear, imonth, iday, ihour, iminute, isecond, itick)
  if (abs(CS%timeunit - 86400.0) < 1.0) then
    reday = REAL(num_days)+ (REAL(start_of_day)/86400.0)
    mesg_intro = "MOM Day "
  else
    reday = REAL(num_days)*(86400.0/CS%timeunit) + &
            REAL(start_of_day)/abs(CS%timeunit)
    mesg_intro = "MOM Time "
  endif
  if (reday < 1.0e8) then ;      write(day_str, '(F12.3)') reday
  elseif (reday < 1.0e11) then ; write(day_str, '(F15.3)') reday
  else ;                         write(day_str, '(ES15.9)') reday ; endif

  if     (n < 1000000)   then ; write(n_str, '(I6)')  n
  elseif (n < 10000000)  then ; write(n_str, '(I7)')  n
  elseif (n < 100000000) then ; write(n_str, '(I8)')  n
  else                        ; write(n_str, '(I10)') n ; endif

  if (date_stamped) then
    write(date_str,'("MOM Date",i7,2("/",i2.2)," ",i2.2,2(":",i2.2))') &
       iyear, imonth, iday, ihour, iminute, isecond
  else
    date_str = trim(mesg_intro)//trim(day_str)
  endif

  if (is_root_pe()) then
    if (CS%use_temperature) then
        write(*,'(A," ",A,": En ",ES12.6, ", MaxCFL ", F8.5, ", Mass ", &
                & ES18.12, ", Salt ", F15.11,", Temp ", F15.11)') &
          trim(date_str), trim(n_str), En_mass, max_CFL(1), mass_tot, salin, temp
    else
        write(*,'(A," ",A,": En ",ES12.6, ", MaxCFL ", F8.5, ", Mass ", &
                & ES18.12)') &
          trim(date_str), trim(n_str), En_mass, max_CFL(1), mass_tot
    endif

    if (CS%use_temperature) then
      write(CS%fileenergy_ascii,'(A,",",A,",", I6,", En ",ES22.16, &
                               &", CFL ", F8.5, ", SL ",&
                               &es11.4,", M ",ES11.5,", S",f8.4,", T",f8.4,&
                               &", Me ",ES9.2,", Se ",ES9.2,", Te ",ES9.2)') &
            trim(n_str), trim(day_str), CS%ntrunc, En_mass, max_CFL(1), &
            -H_0APE(1), mass_tot, salin, temp, mass_anom/mass_tot, salin_anom, &
            temp_anom
    else
      write(CS%fileenergy_ascii,'(A,",",A,",", I6,", En ",ES22.16, &
                                &", CFL ", F8.5, ", SL ",&
                                  &ES11.4,", Mass ",ES11.5,", Me ",ES9.2)') &
            trim(n_str), trim(day_str), CS%ntrunc, En_mass, max_CFL(1), &
            -H_0APE(1), mass_tot, mass_anom/mass_tot
    endif

    if (CS%ntrunc > 0) then
      write(*,'(A," Energy/Mass:",ES12.5," Truncations ",I0)') &
        trim(date_str), En_mass, CS%ntrunc
    endif

    if (CS%write_stocks) then
      write(*,'("    Total Energy: ",Z16.16,ES24.16)') toten, toten
      write(*,'("    Total Mass: ",ES24.16,", Change: ",ES24.16," Error: ",ES12.5," (",ES8.1,")")') &
            mass_tot, mass_chg, mass_anom, mass_anom/mass_tot
      if (CS%use_temperature) then
        if (Salt == 0.) then
          write(*,'("    Total Salt: ",ES24.16,", Change: ",ES24.16," Error: ",ES12.5)') &
              Salt*0.001, Salt_chg*0.001, Salt_anom*0.001
        else
          write(*,'("    Total Salt: ",ES24.16,", Change: ",ES24.16," Error: ",ES12.5," (",ES8.1,")")') &
              Salt*0.001, Salt_chg*0.001, Salt_anom*0.001, Salt_anom/Salt
        endif
        if (Heat == 0.) then
          write(*,'("    Total Heat: ",ES24.16,", Change: ",ES24.16," Error: ",ES12.5)') &
              Heat, Heat_chg, Heat_anom
        else
          write(*,'("    Total Heat: ",ES24.16,", Change: ",ES24.16," Error: ",ES12.5," (",ES8.1,")")') &
              Heat, Heat_chg, Heat_anom, Heat_anom/Heat
        endif
      endif
      do m=1,nTr_stocks

         write(*,'("      Total ",a,": ",ES24.16,X,a)') &
              trim(Tr_names(m)), Tr_stocks(m), trim(Tr_units(m))

         if (Tr_minmax_got(m)) then
           write(*,'(64X,"Global Min:",ES24.16,X,"at: (", f7.2,","f7.2,","f8.2,")"  )') &
                Tr_min(m),Tr_min_x(m),Tr_min_y(m),Tr_min_z(m)
           write(*,'(64X,"Global Max:",ES24.16,X,"at: (", f7.2,","f7.2,","f8.2,")"  )') &
                Tr_max(m),Tr_max_x(m),Tr_max_y(m),Tr_max_z(m)
        endif

      enddo
    endif
  endif

  var = real(CS%ntrunc)
  call write_field(CS%fileenergy_nc, CS%fields(1), var, reday)
  call write_field(CS%fileenergy_nc, CS%fields(2), toten, reday)
  call write_field(CS%fileenergy_nc, CS%fields(3), PE, reday)
  call write_field(CS%fileenergy_nc, CS%fields(4), KE, reday)
  call write_field(CS%fileenergy_nc, CS%fields(5), H_0APE, reday)
  call write_field(CS%fileenergy_nc, CS%fields(6), mass_lay, reday)

  call write_field(CS%fileenergy_nc, CS%fields(7), mass_tot, reday)
  call write_field(CS%fileenergy_nc, CS%fields(8), mass_chg, reday)
  call write_field(CS%fileenergy_nc, CS%fields(9), mass_anom, reday)
  call write_field(CS%fileenergy_nc, CS%fields(10), max_CFL(1), reday)
  call write_field(CS%fileenergy_nc, CS%fields(11), max_CFL(1), reday)
  if (CS%use_temperature) then
    call write_field(CS%fileenergy_nc, CS%fields(12), 0.001*Salt, reday)
    call write_field(CS%fileenergy_nc, CS%fields(13), 0.001*salt_chg, reday)
    call write_field(CS%fileenergy_nc, CS%fields(14), 0.001*salt_anom, reday)
    call write_field(CS%fileenergy_nc, CS%fields(15), Heat, reday)
    call write_field(CS%fileenergy_nc, CS%fields(16), heat_chg, reday)
    call write_field(CS%fileenergy_nc, CS%fields(17), heat_anom, reday)
    do m=1,nTr_stocks
      call write_field(CS%fileenergy_nc, CS%fields(17+m), Tr_stocks(m), reday)
    enddo
  else
    do m=1,nTr_stocks
      call write_field(CS%fileenergy_nc, CS%fields(11+m), Tr_stocks(m), reday)
    enddo
  endif

  call flush_file(CS%fileenergy_nc)

  ! The second (impossible-looking) test looks for a NaN in En_mass.
  if ((En_mass>CS%max_Energy) .or. &
     ((En_mass>CS%max_Energy) .and. (En_mass<CS%max_Energy))) then
    write(mesg,'("Energy per unit mass of ",ES11.4," exceeds ",ES11.4)') &
                  En_mass, CS%max_Energy
    call MOM_error(FATAL, &
      "write_energy : Excessive energy per unit mass or NaNs forced model termination.")
  endif
  if (CS%ntrunc>CS%maxtrunc) then
    call MOM_error(FATAL, "write_energy : Ocean velocity has been truncated too many times.")
  endif
  CS%ntrunc = 0
  CS%previous_calls = CS%previous_calls + 1
  CS%mass_prev = mass_tot ; CS%fresh_water_input = 0.0
  if (CS%use_temperature) then
    CS%salt_prev = Salt ; CS%net_salt_input = 0.0
    CS%heat_prev = Heat ; CS%net_heat_input = 0.0
  endif

  CS%mass_prev_EFP = mass_EFP ; CS%fresh_water_in_EFP = real_to_EFP(0.0)
  if (CS%use_temperature) then
    CS%salt_prev_EFP = Salt_EFP ; CS%net_salt_in_EFP = real_to_EFP(0.0)
    CS%heat_prev_EFP = Heat_EFP ; CS%net_heat_in_EFP = real_to_EFP(0.0)
  endif
end subroutine write_energy

!> This subroutine accumates the net input of volume, salt and heat, through
!! the ocean surface for use in diagnosing conservation.
subroutine accumulate_net_input(fluxes, sfc_state, tv, dt, G, US, CS)
  type(forcing),         intent(in) :: fluxes !< A structure containing pointers to any possible
                                              !! forcing fields.  Unused fields are unallocated.
  type(surface),         intent(in) :: sfc_state !< A structure containing fields that
                                              !! describe the surface state of the ocean.
  type(thermo_var_ptrs), intent(in) :: tv     !< A structure pointing to various
                                              !! thermodynamic variables.
  real,                  intent(in) :: dt     !< The amount of time over which to average [T ~> s].
  type(ocean_grid_type), intent(in) :: G      !< The ocean's grid structure.
  type(unit_scale_type), intent(in) :: US     !< A dimensional unit scaling type
  type(Sum_output_CS),   pointer    :: CS     !< The control structure returned by a previous call
                                              !! to MOM_sum_output_init.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    FW_in, &   ! The net fresh water input, integrated over a timestep [kg].
    salt_in, & ! The total salt added by surface fluxes, integrated
               ! over a time step [ppt kg].
    heat_in    ! The total heat added by surface fluxes, integrated
               ! over a time step [J].
  real :: FW_input   ! The net fresh water input, integrated over a timestep
                     ! and summed over space [kg].
  real :: salt_input ! The total salt added by surface fluxes, integrated
                     ! over a time step and summed over space [ppt kg].
  real :: heat_input ! The total heat added by boundary fluxes, integrated
                     ! over a time step and summed over space [J].
  real :: C_p        ! The heat capacity of seawater [J degC-1 kg-1].
  real :: RZL2_to_kg ! A combination of scaling factors for mass [kg R-1 Z-1 L-2 ~> 1]

  type(EFP_type) :: &
    FW_in_EFP,   & ! Extended fixed point version of FW_input [kg]
    salt_in_EFP, & ! Extended fixed point version of salt_input [ppt kg]
    heat_in_EFP    ! Extended fixed point version of heat_input [J]

  real :: inputs(3)   ! A mixed array for combining the sums
  integer :: i, j, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  C_p = fluxes%C_p
  RZL2_to_kg = US%L_to_m**2*US%R_to_kg_m3*US%Z_to_m

  FW_in(:,:) = 0.0 ; FW_input = 0.0
  if (associated(fluxes%evap)) then
    if (associated(fluxes%lprec) .and. associated(fluxes%fprec)) then
      do j=js,je ; do i=is,ie
        FW_in(i,j) = RZL2_to_kg * dt*G%areaT(i,j)*(fluxes%evap(i,j) + &
            (((fluxes%lprec(i,j) + fluxes%vprec(i,j)) + fluxes%lrunoff(i,j)) + &
              (fluxes%fprec(i,j) + fluxes%frunoff(i,j))))
      enddo ; enddo
    else
      call MOM_error(WARNING, &
        "accumulate_net_input called with associated evap field, but no precip field.")
    endif
  endif

  if (associated(fluxes%seaice_melt)) then ; do j=js,je ; do i=is,ie
    FW_in(i,j) = FW_in(i,j) + RZL2_to_kg*dt * &
                 G%areaT(i,j) * fluxes%seaice_melt(i,j)
  enddo ; enddo ; endif

  salt_in(:,:) = 0.0 ; heat_in(:,:) = 0.0
  if (CS%use_temperature) then

    if (associated(fluxes%sw)) then ; do j=js,je ; do i=is,ie
      heat_in(i,j) = heat_in(i,j) + dt*US%T_to_s*US%L_to_m**2*G%areaT(i,j) * (fluxes%sw(i,j) + &
             (fluxes%lw(i,j) + (fluxes%latent(i,j) + fluxes%sens(i,j))))
    enddo ; enddo ; endif

    if (associated(fluxes%seaice_melt_heat)) then ; do j=js,je ; do i=is,ie
       heat_in(i,j) = heat_in(i,j) + dt*US%T_to_s*US%L_to_m**2*G%areaT(i,j) * fluxes%seaice_melt_heat(i,j)
    enddo ; enddo ; endif

    ! smg: new code
    ! include heat content from water transport across ocean surface
!    if (associated(fluxes%heat_content_lprec)) then ; do j=js,je ; do i=is,ie
!      heat_in(i,j) = heat_in(i,j) + dt*RZL2_to_kg*G%areaT(i,j) * &
!         (fluxes%heat_content_lprec(i,j)   + (fluxes%heat_content_fprec(i,j)   &
!       + (fluxes%heat_content_lrunoff(i,j) + (fluxes%heat_content_frunoff(i,j) &
!       + (fluxes%heat_content_cond(i,j)    + (fluxes%heat_content_vprec(i,j)   &
!       +  fluxes%heat_content_massout(i,j)))))))
!    enddo ; enddo ; endif

    ! smg: old code
    if (associated(tv%TempxPmE)) then
      do j=js,je ; do i=is,ie
        heat_in(i,j) = heat_in(i,j) + (C_p * RZL2_to_kg*G%areaT(i,j)) * tv%TempxPmE(i,j)
      enddo ; enddo
    elseif (associated(fluxes%evap)) then
      do j=js,je ; do i=is,ie
        heat_in(i,j) = heat_in(i,j) + (C_p * sfc_state%SST(i,j)) * FW_in(i,j)
      enddo ; enddo
    endif


    ! The following heat sources may or may not be used.
    if (associated(tv%internal_heat)) then
      do j=js,je ; do i=is,ie
        heat_in(i,j) = heat_in(i,j) + (C_p * US%L_to_m**2*G%areaT(i,j)) * &
                       tv%internal_heat(i,j)
      enddo ; enddo
    endif
    if (associated(tv%frazil)) then ; do j=js,je ; do i=is,ie
      heat_in(i,j) = heat_in(i,j) + US%L_to_m**2*G%areaT(i,j) * tv%frazil(i,j)
    enddo ; enddo ; endif
    if (associated(fluxes%heat_added)) then ; do j=js,je ; do i=is,ie
      heat_in(i,j) = heat_in(i,j) + dt*US%T_to_s*US%L_to_m**2*G%areaT(i,j)*fluxes%heat_added(i,j)
    enddo ; enddo ; endif
!    if (associated(sfc_state%sw_lost)) then ; do j=js,je ; do i=is,ie
!      heat_in(i,j) = heat_in(i,j) - US%L_to_m**2*G%areaT(i,j) * sfc_state%sw_lost(i,j)
!    enddo ; enddo ; endif

    if (associated(fluxes%salt_flux)) then ; do j=js,je ; do i=is,ie
      ! convert salt_flux from kg (salt)/(m^2 s) to ppt * [m s-1].
      salt_in(i,j) = RZL2_to_kg * dt * &
                     G%areaT(i,j)*(1000.0*fluxes%salt_flux(i,j))
    enddo ; enddo ; endif
  endif

  if ((CS%use_temperature) .or. associated(fluxes%lprec) .or. &
      associated(fluxes%evap)) then
    FW_input   = reproducing_sum(FW_in,   EFP_sum=FW_in_EFP)
    heat_input = reproducing_sum(heat_in, EFP_sum=heat_in_EFP)
    salt_input = reproducing_sum(salt_in, EFP_sum=salt_in_EFP)

    CS%fresh_water_input = CS%fresh_water_input + FW_input
    CS%net_salt_input    = CS%net_salt_input    + salt_input
    CS%net_heat_input    = CS%net_heat_input    + heat_input

    CS%fresh_water_in_EFP = CS%fresh_water_in_EFP + FW_in_EFP
    CS%net_salt_in_EFP    = CS%net_salt_in_EFP    + salt_in_EFP
    CS%net_heat_in_EFP    = CS%net_heat_in_EFP    + heat_in_EFP
  endif

end subroutine accumulate_net_input

!>  This subroutine sets up an ordered list of depths, along with the
!! cross sectional areas at each depth and the volume of fluid deeper
!! than each depth.  This might be read from a previously created file
!! or it might be created anew.  (For now only new creation occurs.
subroutine depth_list_setup(G, US, CS)
  type(ocean_grid_type), intent(in) :: G   !< The ocean's grid structure
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  type(Sum_output_CS),   pointer    :: CS  !< The control structure returned by a
                                           !! previous call to MOM_sum_output_init.
  ! Local variables
  integer :: k

  if (CS%read_depth_list) then
    if (file_exists(CS%depth_list_file)) then
      call read_depth_list(G, US, CS, CS%depth_list_file)
    else
      if (is_root_pe()) call MOM_error(WARNING, "depth_list_setup: "// &
        trim(CS%depth_list_file)//" does not exist.  Creating a new file.")
      call create_depth_list(G, CS)

      call write_depth_list(G, US, CS, CS%depth_list_file, CS%list_size+1)
    endif
  else
    call create_depth_list(G, CS)
  endif

  do k=1,G%ke
    CS%lH(k) = CS%list_size
  enddo

end subroutine depth_list_setup

!>  create_depth_list makes an ordered list of depths, along with the cross
!! sectional areas at each depth and the volume of fluid deeper than each depth.
subroutine create_depth_list(G, CS)
  type(ocean_grid_type), intent(in) :: G  !< The ocean's grid structure.
  type(Sum_output_CS),   pointer    :: CS !< The control structure set up in MOM_sum_output_init,
                                          !! in which the ordered depth list is stored.
  ! Local variables
  real, dimension(G%Domain%niglobal*G%Domain%njglobal + 1) :: &
    Dlist, &  !< The global list of bottom depths [Z ~> m].
    AreaList  !< The global list of cell areas [m2].
  integer, dimension(G%Domain%niglobal*G%Domain%njglobal+1) :: &
    indx2     !< The position of an element in the original unsorted list.
  real    :: Dnow  !< The depth now being considered for sorting [Z ~> m].
  real    :: Dprev !< The most recent depth that was considered [Z ~> m].
  real    :: vol   !< The running sum of open volume below a deptn [Z m2 ~> m3].
  real    :: area  !< The open area at the current depth [m2].
  real    :: D_list_prev !< The most recent depth added to the list [Z ~> m].
  logical :: add_to_list !< This depth should be included as an entry on the list.

  integer :: ir, indxt
  integer :: mls, list_size
  integer :: list_pos, i_global, j_global
  integer :: i, j, k, kl

  mls = G%Domain%niglobal*G%Domain%njglobal

! Need to collect the global data from compute domains to a 1D array for sorting.
  Dlist(:) = 0.0
  Arealist(:) = 0.0
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    ! Set global indices that start the global domain at 1 (Fortran convention).
    j_global = j + G%jdg_offset - (G%jsg-1)
    i_global = i + G%idg_offset - (G%isg-1)

    list_pos = (j_global-1)*G%Domain%niglobal + i_global
    Dlist(list_pos) = G%bathyT(i,j)
    Arealist(list_pos) = G%mask2dT(i,j) * G%US%L_to_m**2*G%areaT(i,j)
  enddo ; enddo

  ! These sums reproduce across PEs because the arrays are only nonzero on one PE.
  call sum_across_PEs(Dlist, mls+1)
  call sum_across_PEs(Arealist, mls+1)

  do j=1,mls+1 ; indx2(j) = j ; enddo
  k = mls / 2  + 1 ; ir = mls
  do
    if (k > 1) then
      k = k - 1
      indxt = indx2(k)
      Dnow = Dlist(indxt)
    else
      indxt = indx2(ir)
      Dnow = Dlist(indxt)
      indx2(ir) = indx2(1)
      ir = ir - 1
      if (ir == 1) then ; indx2(1) = indxt ; exit ; endif
    endif
    i=k ; j=k*2
    do ; if (j > ir) exit
      if (j < ir .AND. Dlist(indx2(j)) < Dlist(indx2(j+1))) j = j + 1
      if (Dnow < Dlist(indx2(j))) then ; indx2(i) = indx2(j) ; i = j ; j = j + i
      else ; j = ir+1 ; endif
    enddo
    indx2(i) = indxt
  enddo

!  At this point, the lists should perhaps be culled to save memory.
! Culling: (1) identical depths (e.g. land) - take the last one.
!          (2) the topmost and bottommost depths are always saved.
!          (3) very close depths
!          (4) equal volume changes.

  ! Count the unique elements in the list.
  D_list_prev = Dlist(indx2(mls))
  list_size = 2
  do k=mls-1,1,-1
    if (Dlist(indx2(k)) < D_list_prev-CS%D_list_min_inc) then
      list_size = list_size + 1
      D_list_prev = Dlist(indx2(k))
    endif
  enddo

  CS%list_size = list_size
  allocate(CS%DL(CS%list_size+1))

  vol = 0.0 ; area = 0.0
  Dprev = Dlist(indx2(mls))
  D_list_prev = Dprev

  kl = 0
  do k=mls,1,-1
    i = indx2(k)
    vol = vol + area * (Dprev - Dlist(i))
    area = area + AreaList(i)

    add_to_list = .false.
    if ((kl == 0) .or. (k==1)) then
      add_to_list = .true.
    elseif (Dlist(indx2(k-1)) < D_list_prev-CS%D_list_min_inc) then
      add_to_list = .true.
      D_list_prev = Dlist(indx2(k-1))
    endif

    if (add_to_list) then
      kl = kl+1
      CS%DL(kl)%depth = Dlist(i)
      CS%DL(kl)%area = area
      CS%DL(kl)%vol_below = vol
    endif
    Dprev = Dlist(i)
  enddo

  do while (kl < list_size)
    ! I don't understand why this is needed... RWH
    kl = kl+1
    CS%DL(kl)%vol_below = CS%DL(kl-1)%vol_below * 1.000001
    CS%DL(kl)%area = CS%DL(kl-1)%area
    CS%DL(kl)%depth = CS%DL(kl-1)%depth
  enddo

  CS%DL(CS%list_size+1)%vol_below = CS%DL(CS%list_size)%vol_below * 1000.0
  CS%DL(CS%list_size+1)%area = CS%DL(CS%list_size)%area
  CS%DL(CS%list_size+1)%depth = CS%DL(CS%list_size)%depth

end subroutine create_depth_list

!> This subroutine writes out the depth list to the specified file.
subroutine write_depth_list(G, US, CS, filename, list_size)
  type(ocean_grid_type), intent(in) :: G   !< The ocean's grid structure.
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  type(Sum_output_CS),   pointer    :: CS  !< The control structure returned by a
                                           !! previous call to MOM_sum_output_init.
  character(len=*),      intent(in) :: filename !< The path to the depth list file to write.
  integer,               intent(in) :: list_size !< The size of the depth list.
  ! Local variables
  real, allocatable :: tmp(:)
  integer :: ncid, dimid(1), Did, Aid, Vid, status, k
  character(len=16) :: depth_chksum, area_chksum

  ! All ranks are required to compute the global checksum
  call get_depth_list_checksums(G, depth_chksum, area_chksum)

  if (.not.is_root_pe()) return

  allocate(tmp(list_size)) ; tmp(:) = 0.0

  status = NF90_CREATE(filename, 0, ncid)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING, filename//trim(NF90_STRERROR(status)))
    return
  endif

  status = NF90_DEF_DIM(ncid, "list", list_size, dimid(1))
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//trim(NF90_STRERROR(status)))

  status = NF90_DEF_VAR(ncid, "depth", NF90_DOUBLE, dimid, Did)
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" depth "//trim(NF90_STRERROR(status)))
  status = NF90_PUT_ATT(ncid, Did, "long_name", "Sorted depth")
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" depth "//trim(NF90_STRERROR(status)))
  status = NF90_PUT_ATT(ncid, Did, "units", "m")
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" depth "//trim(NF90_STRERROR(status)))

  status = NF90_DEF_VAR(ncid, "area", NF90_DOUBLE, dimid, Aid)
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" area "//trim(NF90_STRERROR(status)))
  status = NF90_PUT_ATT(ncid, Aid, "long_name", "Open area at depth")
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" area "//trim(NF90_STRERROR(status)))
  status = NF90_PUT_ATT(ncid, Aid, "units", "m2")
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" area "//trim(NF90_STRERROR(status)))

  status = NF90_DEF_VAR(ncid, "vol_below", NF90_DOUBLE, dimid, Vid)
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" vol_below "//trim(NF90_STRERROR(status)))
  status = NF90_PUT_ATT(ncid, Vid, "long_name", "Open volume below depth")
   if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" vol_below "//trim(NF90_STRERROR(status)))
  status = NF90_PUT_ATT(ncid, Vid, "units", "m3")
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" vol_below "//trim(NF90_STRERROR(status)))

  ! Dependency checksums
  status = NF90_PUT_ATT(ncid, NF90_GLOBAL, depth_chksum_attr, depth_chksum)
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" "//depth_chksum_attr//" "//trim(NF90_STRERROR(status)))

  status = NF90_PUT_ATT(ncid, NF90_GLOBAL, area_chksum_attr, area_chksum)
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" "//area_chksum_attr//" "//trim(NF90_STRERROR(status)))

  status = NF90_ENDDEF(ncid)
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//trim(NF90_STRERROR(status)))

  do k=1,list_size ; tmp(k) = US%Z_to_m*CS%DL(k)%depth ; enddo
  status = NF90_PUT_VAR(ncid, Did, tmp)
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" depth "//trim(NF90_STRERROR(status)))

  do k=1,list_size ; tmp(k) = CS%DL(k)%area ; enddo
  status = NF90_PUT_VAR(ncid, Aid, tmp)
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" area "//trim(NF90_STRERROR(status)))

  do k=1,list_size ; tmp(k) = US%Z_to_m*CS%DL(k)%vol_below ; enddo
  status = NF90_PUT_VAR(ncid, Vid, tmp)
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//" vol_below "//trim(NF90_STRERROR(status)))

  status = NF90_CLOSE(ncid)
  if (status /= NF90_NOERR) call MOM_error(WARNING, &
      filename//trim(NF90_STRERROR(status)))

end subroutine write_depth_list

!> This subroutine reads in the depth list to the specified file
!! and allocates and sets up CS%DL and CS%list_size .
subroutine read_depth_list(G, US, CS, filename)
  type(ocean_grid_type), intent(in) :: G   !< The ocean's grid structure
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  type(Sum_output_CS),   pointer    :: CS  !< The control structure returned by a
                                           !! previous call to MOM_sum_output_init.
  character(len=*),      intent(in) :: filename !< The path to the depth list file to read.
  ! Local variables
  character(len=32) :: mdl
  character(len=240) :: var_name, var_msg
  real, allocatable :: tmp(:)
  integer :: ncid, status, varid, list_size, k
  integer :: ndim, len, var_dim_ids(NF90_MAX_VAR_DIMS)
  character(len=16) :: depth_file_chksum, depth_grid_chksum
  character(len=16) :: area_file_chksum, area_grid_chksum
  integer :: depth_attr_status, area_attr_status

  mdl = "MOM_sum_output read_depth_list:"

  status = NF90_OPEN(filename, NF90_NOWRITE, ncid)
  if (status /= NF90_NOERR) then
    call MOM_error(FATAL,mdl//" Difficulties opening "//trim(filename)// &
        " - "//trim(NF90_STRERROR(status)))
  endif

  ! Check bathymetric consistency
  depth_attr_status = NF90_GET_ATT(ncid, NF90_GLOBAL, depth_chksum_attr, &
                                   depth_file_chksum)
  area_attr_status = NF90_GET_ATT(ncid, NF90_GLOBAL, area_chksum_attr, &
                                  area_file_chksum)

  if (any([depth_attr_status, area_attr_status] == NF90_ENOTATT)) then
    var_msg = trim(CS%depth_list_file) // " checksums are missing;"
    if (CS%require_depth_list_chksum) then
      call MOM_error(FATAL, trim(var_msg) // " aborting.")
    elseif (CS%update_depth_list_chksum) then
      call MOM_error(WARNING, trim(var_msg) // " updating file.")
      call create_depth_list(G, CS)
      call write_depth_list(G, US, CS, CS%depth_list_file, CS%list_size+1)
      return
    else
      call MOM_error(WARNING, &
        trim(var_msg) // " some diagnostics may not be reproducible.")
    endif
  else
    ! Validate netCDF call
    if (depth_attr_status /= NF90_NOERR) then
      var_msg = mdl // "Failed to read " // trim(filename) // ":" &
                // depth_chksum_attr
      call MOM_error(FATAL, &
        trim(var_msg) // " - " // NF90_STRERROR(depth_attr_status))
    endif

    if (area_attr_status /= NF90_NOERR) then
      var_msg = mdl // "Failed to read " // trim(filename) // ":" &
                // area_chksum_attr
      call MOM_error(FATAL, &
        trim(var_msg) // " - " // NF90_STRERROR(area_attr_status))
    endif

    call get_depth_list_checksums(G, depth_grid_chksum, area_grid_chksum)

    if (depth_grid_chksum /= depth_file_chksum &
            .or. area_grid_chksum /= area_file_chksum) then
      var_msg = trim(CS%depth_list_file) // " checksums do not match;"
      if (CS%require_depth_list_chksum) then
        call MOM_error(FATAL, trim(var_msg) // " aborting.")
      elseif (CS%update_depth_list_chksum) then
        call MOM_error(WARNING, trim(var_msg) // " updating file.")
        call create_depth_list(G, CS)
        call write_depth_list(G, US, CS, CS%depth_list_file, CS%list_size+1)
        return
      else
        call MOM_error(WARNING, &
          trim(var_msg) // " some diagnostics may not be reproducible.")
      endif
    endif
  endif

  var_name = "depth"
  var_msg = trim(var_name)//" in "//trim(filename)//" - "
  status = NF90_INQ_VARID(ncid, var_name, varid)
  if (status /= NF90_NOERR) call MOM_error(FATAL,mdl// &
        " Difficulties finding variable "//trim(var_msg)//&
        trim(NF90_STRERROR(status)))

  status = NF90_INQUIRE_VARIABLE(ncid, varid, ndims=ndim, dimids=var_dim_ids)
  if (status /= NF90_NOERR) then
    call MOM_ERROR(FATAL,mdl//" cannot inquire about "//trim(var_msg)//&
        trim(NF90_STRERROR(status)))
  elseif (ndim > 1) then
    call MOM_ERROR(FATAL,mdl//" "//trim(var_msg)//&
         " has too many or too few dimensions.")
  endif

  ! Get the length of the list.
  status = NF90_INQUIRE_DIMENSION(ncid, var_dim_ids(1), len=list_size)
  if (status /= NF90_NOERR) call MOM_ERROR(FATAL,mdl// &
        " cannot inquire about dimension(1) of "//trim(var_msg)//&
        trim(NF90_STRERROR(status)))

  CS%list_size = list_size-1
  allocate(CS%DL(list_size))
  allocate(tmp(list_size))

  status = NF90_GET_VAR(ncid, varid, tmp)
  if (status /= NF90_NOERR) call MOM_error(FATAL,mdl// &
        " Difficulties reading variable "//trim(var_msg)//&
        trim(NF90_STRERROR(status)))

  do k=1,list_size ; CS%DL(k)%depth = US%m_to_Z*tmp(k) ; enddo

  var_name = "area"
  var_msg = trim(var_name)//" in "//trim(filename)//" - "
  status = NF90_INQ_VARID(ncid, var_name, varid)
  if (status /= NF90_NOERR) call MOM_error(FATAL,mdl// &
        " Difficulties finding variable "//trim(var_msg)//&
        trim(NF90_STRERROR(status)))
  status = NF90_GET_VAR(ncid, varid, tmp)
  if (status /= NF90_NOERR) call MOM_error(FATAL,mdl// &
        " Difficulties reading variable "//trim(var_msg)//&
        trim(NF90_STRERROR(status)))

  do k=1,list_size ; CS%DL(k)%area = tmp(k) ; enddo

  var_name = "vol_below"
  var_msg = trim(var_name)//" in "//trim(filename)
  status = NF90_INQ_VARID(ncid, var_name, varid)
  if (status /= NF90_NOERR) call MOM_error(FATAL,mdl// &
        " Difficulties finding variable "//trim(var_msg)//&
        trim(NF90_STRERROR(status)))
  status = NF90_GET_VAR(ncid, varid, tmp)
  if (status /= NF90_NOERR) call MOM_error(FATAL,mdl// &
        " Difficulties reading variable "//trim(var_msg)//&
        trim(NF90_STRERROR(status)))

  do k=1,list_size ; CS%DL(k)%vol_below = US%m_to_Z*tmp(k) ; enddo

  status = NF90_CLOSE(ncid)
  if (status /= NF90_NOERR) call MOM_error(WARNING, mdl// &
    " Difficulties closing "//trim(filename)//" - "//trim(NF90_STRERROR(status)))

  deallocate(tmp)

end subroutine read_depth_list


!> Return the checksums required to verify DEPTH_LIST_FILE contents.
!!
!! This function computes checksums for the bathymetry (G%bathyT) and masked
!! area (mask2dT * areaT) fields of the model grid G, which are used to compute
!! the depth list.  A difference in checksum indicates that a different method
!! was used to compute the grid data, and that any results using the depth
!! list, such as APE, will not be reproducible.
!!
!! Checksums are saved as hexadecimal strings, in order to avoid potential
!! datatype issues with netCDF attributes.
subroutine get_depth_list_checksums(G, depth_chksum, area_chksum)
  type(ocean_grid_type), intent(in) :: G          !< Ocean grid structure
  character(len=16), intent(out) :: depth_chksum  !< Depth checksum hexstring
  character(len=16), intent(out) :: area_chksum   !< Area checksum hexstring

  integer :: i, j
  real, allocatable :: field(:,:)

  allocate(field(G%isc:G%iec, G%jsc:G%jec))

  ! Depth checksum
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    field(i,j) = G%bathyT(i,j)
  enddo ; enddo
  write(depth_chksum, '(Z16)') mpp_chksum(field(:,:))

  ! Area checksum
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    field(i,j) = G%mask2dT(i,j) * G%US%L_to_m**2*G%areaT(i,j)
  enddo ; enddo
  write(area_chksum, '(Z16)') mpp_chksum(field(:,:))

  deallocate(field)
end subroutine get_depth_list_checksums

!> \namespace mom_sum_output
!!
!! By Robert Hallberg, April 1994 - June 2002
!!
!!   This file contains the subroutine (write_energy) that writes
!! horizontally integrated quantities, such as energies and layer
!! volumes, and other summary information to an output file.  Some
!! of these quantities (APE or resting interface height) are defined
!! relative to the global histogram of topography.  The subroutine
!! that compiles that histogram (depth_list_setup) is also included
!! in this file.
!!
!!   In addition, if the number of velocity truncations since the
!! previous call to write_energy exceeds maxtrunc or the total energy
!! exceeds a very large threshold, a fatal termination is triggered.

end module MOM_sum_output
