module MOM_OCMIP2_CFC
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
!*  By Robert Hallberg, 2007                                           *
!*                                                                     *
!*    This file contains an example of the code that is needed to set  *
!*  up and use CFC-11 and CFC-12 in a fully coupled or ice-ocean model *
!*  context. There are 5 subroutines in this file.                     *
!*                                                                     *
!*    register_OCMIP2_CFC determines if the module is going to work,   *
!*  then makes several calls registering tracers to be advected and    *
!*  read from a restart file. it also sets various run-time parameters *
!*  for this module and sets up a "control structure" (CS) to store    *
!*  all information for this module.                                   *
!*                                                                     *
!*    initialize_OCMIP2_CFC initializes this modules arrays if they    *
!*  have not been found in a restart file.  It also determines which   *
!*  diagnostics will need to be calculated.                            *
!*                                                                     *
!*    OCMIP2_CFC_column_physics updates the CFC concentrations,        *
!*  applying everthing but horizontal advection and diffusion.         *
!*  Surface fluxes are applied inside an implicit vertical advection   *
!*  and diffusion tridiagonal solver, and any interior sources and     *
!*  sinks (not applicable for CFCs) would also be applied here.  This  *
!*  subroutine also sends out any requested interior diagnostics.      *
!*                                                                     *
!*    OCMIP2_CFC_surface_state calculates the information required     *
!*  from the ocean for the FMS coupler to calculate CFC fluxes.        *
!*                                                                     *
!*    OCMIP2_CFC_end deallocates the persistent run-time memory used   *
!*  by this module.                                                    *
!*                                                                     *
!*     A small fragment of the horizontal grid is shown below:         *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, tr_ady, tr_dfy                        *
!*    j    x ^ x ^ x   At >:  u, tr_adx, tr_dfx                        *
!*    j    > o > o >   At o:  h, tr, CFC11, CFC12                      *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl
use MOM_diag_to_Z, only : register_Z_tracer, diag_to_Z_CS
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_hor_index, only : hor_index_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : register_restart_field, query_initialized, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type, get_time
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_registry, only : add_tracer_diagnostics, add_tracer_OBC_values
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init, only : tracer_Z_init
use MOM_variables, only : surface
use MOM_verticalGrid, only : verticalGrid_type

use coupler_util, only : extract_coupler_values, set_coupler_values
use coupler_util, only : ind_flux, ind_alpha, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

public register_OCMIP2_CFC, initialize_OCMIP2_CFC
public OCMIP2_CFC_column_physics, OCMIP2_CFC_surface_state
public OCMIP2_CFC_stock, OCMIP2_CFC_end


! NTR is the number of tracers in this module.
integer, parameter :: NTR = 2

type p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d

type, public :: OCMIP2_CFC_CS ; private
  character(len=200) :: IC_file ! The file in which the CFC initial values can
                    ! be found, or an empty string for internal initilaization.
  logical :: Z_IC_file ! If true, the IC_file is in Z-space.  The default is false..
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL()
  real, pointer, dimension(:,:,:) :: &
    CFC11 => NULL(), &     ! The CFC11 concentration in mol m-3.
    CFC12 => NULL(), &     ! The CFC12 concentration in mol m-3.
    CFC11_aux => NULL(), & ! The CFC11 and CFC12 concentrations, in mol m-3,
    CFC12_aux => NULL()    ! with values of thin layers masked out.
  ! In the following variables a suffix of _11 refers to CFC11 and _12 to CFC12.
  real :: a1_11, a2_11, a3_11, a4_11   ! Coefficients in the calculation of the
  real :: a1_12, a2_12, a3_12, a4_12   ! CFC11 and CFC12 Schmidt numbers, in
                                       ! units of ND, degC-1, degC-2, degC-3.
  real :: d1_11, d2_11, d3_11, d4_11   ! Coefficients in the calculation of the
  real :: d1_12, d2_12, d3_12, d4_12   ! CFC11 and CFC12 solubilities, in units
                                       ! of ND, K-1, log(K)^-1, K-2.
  real :: e1_11, e2_11, e3_11          ! More coefficients in the calculation of
  real :: e1_12, e2_12, e3_12          ! the CFC11 and CFC12 solubilities, in
                                       ! units of PSU-1, PSU-1 K-1, PSU-1 K-2.
  type(p3d), dimension(NTR) :: &
    tr_adx, &       ! Tracer zonal advective fluxes in mol s-1.
    tr_ady, &       ! Tracer meridional advective fluxes in mol s-1.
    tr_dfx, &       ! Tracer zonal diffusive fluxes in mol s-1.
    tr_dfy          ! Tracer meridional diffusive fluxes in mol s-1.
  real :: CFC11_IC_val = 0.0    ! The initial value assigned to CFC11.
  real :: CFC12_IC_val = 0.0    ! The initial value assigned to CFC12.
  real :: CFC11_land_val = -1.0 ! The values of CFC11 and CFC12 used where
  real :: CFC12_land_val = -1.0 ! land is masked out.
  logical :: mask_tracers  ! If true, tracers are masked out in massless layers.
  logical :: tracers_may_reinit  ! If true, tracers may go through the
                           ! initialization code if they are not found in the
                           ! restart files.
  character(len=16) :: CFC11_name, CFC12_name ! Variable names.

  integer :: ind_cfc_11_flux  ! Indices returned by aof_set_coupler_flux that
  integer :: ind_cfc_12_flux  ! are used to pack and unpack surface boundary
                              ! condition arrays.

  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()
  integer :: id_CFC11, id_CFC12
  integer, dimension(NTR) :: id_tr_adx = -1, id_tr_ady = -1
  integer, dimension(NTR) :: id_tr_dfx = -1, id_tr_dfy = -1

  ! The following vardesc types contain a package of metadata about each tracer.
  type(vardesc) :: CFC11_desc, CFC12_desc
end type OCMIP2_CFC_CS

contains

function register_OCMIP2_CFC(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),    intent(in) :: HI
  type(verticalGrid_type), intent(in) :: GV
  type(param_file_type),   intent(in) :: param_file
  type(OCMIP2_CFC_CS),     pointer    :: CS
  type(tracer_registry_type), pointer :: tr_Reg
  type(MOM_restart_CS),    pointer    :: restart_CS
! This subroutine is used to register tracer fields and subroutines
! to be used with MOM.
! Arguments: HI - A horizontal index type structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in/out)  tr_Reg - A pointer to the tracer registry.
!  (in)      restart_CS - A pointer to the restart control structure.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_OCMIP2_CFC" ! This module's name.
  character(len=200) :: inputdir ! The directory where NetCDF input files are.
  ! These can be overridden later in via the field manager?
  character(len=128) :: default_ice_restart_file = 'ice_ocmip2_cfc.res.nc'
  character(len=128) :: default_ocean_restart_file = 'ocmip2_cfc.res.nc'
  real, dimension(:,:,:), pointer :: tr_ptr
  real :: a11_dflt(4), a12_dflt(4) ! Default values of the various coefficients
  real :: d11_dflt(4), d12_dflt(4) ! In the expressions for the solubility and
  real :: e11_dflt(3), e12_dflt(3) ! Schmidt numbers.
  logical :: register_OCMIP2_CFC
  integer :: isd, ied, jsd, jed, nz, m

  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "register_OCMIP2_CFC called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! These calls obtain the indices for the CFC11 and CFC12 flux coupling.
  CS%ind_cfc_11_flux = aof_set_coupler_flux('cfc_11_flux', &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2', &
       param = (/ 9.36e-07, 9.7561e-06 /), &
       ice_restart_file = default_ice_restart_file, &
       ocean_restart_file = default_ocean_restart_file, &
       caller = "register_OCMIP2_CFC")
  CS%ind_cfc_12_flux = aof_set_coupler_flux('cfc_12_flux', &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2', &
       param = (/ 9.36e-07, 9.7561e-06 /), &
       ice_restart_file = default_ice_restart_file, &
       ocean_restart_file = default_ocean_restart_file, &
       caller = "register_OCMIP2_CFC")
  if ((CS%ind_cfc_11_flux < 0) .or. (CS%ind_cfc_11_flux < 0)) then
    ! This is most likely to happen with the dummy version of aof_set_coupler_flux
    ! used in ocean-only runs.
    call MOM_ERROR(WARNING, "CFCs are currently only set up to be run in " // &
                   " coupled model configurations, and will be disabled.")
    deallocate(CS)
    register_OCMIP2_CFC = .false.
    return
  endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "CFC_IC_FILE", CS%IC_file, &
                 "The file in which the CFC initial values can be \n"//&
                 "found, or an empty string for internal initialization.", &
                 default=" ")
  if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file,'/') == 0)) then
    ! Add the directory if CS%IC_file is not already a complete path.
    call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    call log_param(param_file, mod, "INPUTDIR/CFC_IC_FILE", CS%IC_file)
  endif
  call get_param(param_file, mod, "CFC_IC_FILE_IS_Z", CS%Z_IC_file, &
                 "If true, CFC_IC_FILE is in depth space, not layer space", &
                 default=.false.)
  call get_param(param_file, mod, "MASK_MASSLESS_TRACERS", CS%mask_tracers, &
                 "If true, the tracers are masked out in massless layer. \n"//&
                 "This can be a problem with time-averages.", default=.false.)
  call get_param(param_file, mod, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code \n"//&
                 "if they are not found in the restart files.  Otherwise \n"//&
                 "it is a fatal error if tracers are not found in the \n"//&
                 "restart files of a restarted run.", default=.false.)

  !   The following vardesc types contain a package of metadata about each tracer,
  ! including, the name; units; longname; and grid information.
  CS%CFC11_name = "CFC11" ; CS%CFC12_name = "CFC12"
  CS%CFC11_desc = var_desc(CS%CFC11_name,"mol m-3","CFC-11 Concentration", caller=mod)
  CS%CFC12_desc = var_desc(CS%CFC12_name,"mol m-3","CFC-12 Concentration", caller=mod)

  allocate(CS%CFC11(isd:ied,jsd:jed,nz)) ; CS%CFC11(:,:,:) = 0.0
  allocate(CS%CFC12(isd:ied,jsd:jed,nz)) ; CS%CFC12(:,:,:) = 0.0
  if (CS%mask_tracers) then
    allocate(CS%CFC11_aux(isd:ied,jsd:jed,nz)) ; CS%CFC11_aux(:,:,:) = 0.0
    allocate(CS%CFC12_aux(isd:ied,jsd:jed,nz)) ; CS%CFC12_aux(:,:,:) = 0.0
  endif

  ! This pointer assignment is needed to force the compiler not to do a copy in
  ! the registration calls.  Curses on the designers and implementers of F90.
  tr_ptr => CS%CFC11
  ! Register CFC11 for the restart file.
  call register_restart_field(tr_ptr, CS%CFC11_desc, &
                              .not.CS%tracers_may_reinit, restart_CS)
  ! Register CFC11 for horizontal advection & diffusion.
  call register_tracer(tr_ptr, CS%CFC11_desc, param_file, HI, GV, tr_Reg, &
                       tr_desc_ptr=CS%CFC11_desc)
  ! Do the same for CFC12
  tr_ptr => CS%CFC12
  call register_restart_field(tr_ptr, CS%CFC12_desc, &
                              .not.CS%tracers_may_reinit, restart_CS)
  call register_tracer(tr_ptr, CS%CFC12_desc, param_file, HI, GV, tr_Reg, &
                       tr_desc_ptr=CS%CFC12_desc)

  ! Set and read the various empirical coefficients.

!-----------------------------------------------------------------------
! Default Schmidt number coefficients for CFC11 (_11) and CFC12 (_12) are given
! by Zheng et al (1998), JGR vol 103, C1.
!-----------------------------------------------------------------------
  a11_dflt(:) = (/ 3501.8, -210.31,  6.1851, -0.07513 /)
  a12_dflt(:) = (/ 3845.4, -228.95,  6.1908, -0.06743 /)
  call get_param(param_file, mod, "CFC11_A1", CS%a1_11, &
                 "A coefficient in the Schmidt number of CFC11.", &
                 units="nondim", default=a11_dflt(1))
  call get_param(param_file, mod, "CFC11_A2", CS%a2_11, &
                 "A coefficient in the Schmidt number of CFC11.", &
                 units="degC-1", default=a11_dflt(2))
  call get_param(param_file, mod, "CFC11_A3", CS%a3_11, &
                 "A coefficient in the Schmidt number of CFC11.", &
                 units="degC-2", default=a11_dflt(3))
  call get_param(param_file, mod, "CFC11_A4", CS%a4_11, &
                 "A coefficient in the Schmidt number of CFC11.", &
                 units="degC-3", default=a11_dflt(4))

  call get_param(param_file, mod, "CFC12_A1", CS%a1_12, &
                 "A coefficient in the Schmidt number of CFC12.", &
                 units="nondim", default=a12_dflt(1))
  call get_param(param_file, mod, "CFC12_A2", CS%a2_12, &
                 "A coefficient in the Schmidt number of CFC12.", &
                 units="degC-1", default=a12_dflt(2))
  call get_param(param_file, mod, "CFC12_A3", CS%a3_12, &
                 "A coefficient in the Schmidt number of CFC12.", &
                 units="degC-2", default=a12_dflt(3))
  call get_param(param_file, mod, "CFC12_A4", CS%a4_12, &
                 "A coefficient in the Schmidt number of CFC12.", &
                 units="degC-3", default=a12_dflt(4))

!-----------------------------------------------------------------------
! Solubility coefficients for alpha in mol/l/atm for CFC11 (_11) and CFC12 (_12)
! after Warner and Weiss (1985) DSR, vol 32.
!-----------------------------------------------------------------------
  d11_dflt(:) = (/ -229.9261, 319.6552, 119.4471, -1.39165 /)
  e11_dflt(:) = (/ -0.142382, 0.091459, -0.0157274 /)
  d12_dflt(:) = (/ -218.0971, 298.9702, 113.8049, -1.39165 /)
  e12_dflt(:) = (/ -0.143566, 0.091015, -0.0153924 /)

  call get_param(param_file, mod, "CFC11_D1", CS%d1_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="none", default=d11_dflt(1))
  call get_param(param_file, mod, "CFC11_D2", CS%d2_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="hK", default=d11_dflt(2))
  call get_param(param_file, mod, "CFC11_D3", CS%d3_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="none", default=d11_dflt(3))
  call get_param(param_file, mod, "CFC11_D4", CS%d4_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="hK-2", default=d11_dflt(4))
  call get_param(param_file, mod, "CFC11_E1", CS%e1_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="PSU-1", default=e11_dflt(1))
  call get_param(param_file, mod, "CFC11_E2", CS%e2_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="PSU-1 hK-1", default=e11_dflt(2))
  call get_param(param_file, mod, "CFC11_E3", CS%e3_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="PSU-1 hK-2", default=e11_dflt(3))

  call get_param(param_file, mod, "CFC12_D1", CS%d1_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="none", default=d12_dflt(1))
  call get_param(param_file, mod, "CFC12_D2", CS%d2_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="hK", default=d12_dflt(2))
  call get_param(param_file, mod, "CFC12_D3", CS%d3_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="none", default=d12_dflt(3))
  call get_param(param_file, mod, "CFC12_D4", CS%d4_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="hK-2", default=d12_dflt(4))
  call get_param(param_file, mod, "CFC12_E1", CS%e1_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="PSU-1", default=e12_dflt(1))
  call get_param(param_file, mod, "CFC12_E2", CS%e2_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="PSU-1 hK-1", default=e12_dflt(2))
  call get_param(param_file, mod, "CFC12_E3", CS%e3_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="PSU-1 hK-2", default=e12_dflt(3))

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS

  register_OCMIP2_CFC = .true.
end function register_OCMIP2_CFC

subroutine initialize_OCMIP2_CFC(restart, day, G, GV, h, diag, OBC, CS, &
                                 sponge_CSp, diag_to_Z_CSp)
  logical,                               intent(in) :: restart
  type(time_type), target,               intent(in) :: day
  type(ocean_grid_type),                 intent(in) :: G
  type(verticalGrid_type),               intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
  type(diag_ctrl), target,               intent(in) :: diag
  type(ocean_OBC_type),                  pointer    :: OBC
  type(OCMIP2_CFC_CS),                   pointer    :: CS
  type(sponge_CS),                       pointer    :: sponge_CSp
  type(diag_to_Z_CS),                    pointer    :: diag_to_Z_CSp
!   This subroutine initializes the NTR tracer fields in tr(:,:,:,:)
! and it sets up the tracer output.

! Arguments: restart - .true. if the fields have already been read from
!                     a restart file.
!  (in)      day - Time of the start of the run.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in)      OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in/out)  CS - The control structure returned by a previous call to
!                 register_OCMIP2_CFC.
!  (in/out)  sponge_CSp - A pointer to the control structure for the sponges, if
!                         they are in use.  Otherwise this may be unassociated.
!  (in/out)  diag_to_Z_Csp - A pointer to the control structure for diagnostics
!                            in depth space.
  logical :: from_file = .false.
  character(len=16) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for tracer fluxes.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB

  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%Time => day
  CS%diag => diag

  if (.not.restart .or. (CS%tracers_may_reinit .and. &
      .not.query_initialized(CS%CFC11, CS%CFC11_name, CS%restart_CSp))) &
    call init_tracer_CFC(h, CS%CFC11, CS%CFC11_name, CS%CFC11_land_val, &
                         CS%CFC11_IC_val, G, CS)

  if (.not.restart .or. (CS%tracers_may_reinit .and. &
      .not.query_initialized(CS%CFC12, CS%CFC12_name, CS%restart_CSp))) &
    call init_tracer_CFC(h, CS%CFC12, CS%CFC12_name, CS%CFC12_land_val, &
                         CS%CFC12_IC_val, G, CS)

  if (associated(OBC)) then
  ! By default, all tracers have 0 concentration in their inflows. This may
  ! make the following calls are unnecessary.
  !  call add_tracer_OBC_values(trim(CS%CFC11_desc%name), CS%tr_Reg, 0.0)
  !  call add_tracer_OBC_values(trim(CS%CFC12_desc%name), CS%tr_Reg, 0.0)
  endif


  ! This needs to be changed if the units of tracer are changed above.
  if (GV%Boussinesq) then ; flux_units = "mol s-1"
  else ; flux_units = "mol m-3 kg s-1" ; endif

  do m=1,NTR
    ! Register the tracer advective and diffusive fluxes for potential
    ! diagnostic output.
    if (m==1) then
      ! Register CFC11 for potential diagnostic output.
      call query_vardesc(CS%CFC11_desc, name, units=units, longname=longname, &
                         caller="initialize_OCMIP2_CFC")
      CS%id_CFC11 = register_diag_field("ocean_model", trim(name), CS%diag%axesTL, &
          day, trim(longname) , trim(units))
      call register_Z_tracer(CS%CFC11, trim(name), longname, units, &
                             day, G, diag_to_Z_CSp)
    elseif (m==2) then
      ! Register CFC12 for potential diagnostic output.
      call query_vardesc(CS%CFC12_desc, name, units=units, longname=longname, &
                         caller="initialize_OCMIP2_CFC")
      CS%id_CFC12 = register_diag_field("ocean_model", trim(name), CS%diag%axesTL, &
          day, trim(longname) , trim(units))
      call register_Z_tracer(CS%CFC12, trim(name), longname, units, &
                             day, G, diag_to_Z_CSp)
    else
      call MOM_error(FATAL,"initialize_OCMIP2_CFC is only set up to work"//&
                           "with NTR <= 2.")
    endif

    CS%id_tr_adx(m) = register_diag_field("ocean_model", trim(name)//"_adx", &
        CS%diag%axesCuL, day, trim(longname)//" advective zonal flux" , &
        trim(flux_units))
    CS%id_tr_ady(m) = register_diag_field("ocean_model", trim(name)//"_ady", &
        CS%diag%axesCvL, day, trim(longname)//" advective meridional flux" , &
        trim(flux_units))
    CS%id_tr_dfx(m) = register_diag_field("ocean_model", trim(name)//"_dfx", &
        CS%diag%axesCuL, day, trim(longname)//" diffusive zonal flux" , &
        trim(flux_units))
    CS%id_tr_dfy(m) = register_diag_field("ocean_model", trim(name)//"_dfy", &
        CS%diag%axesCvL, day, trim(longname)//" diffusive zonal flux" , &
        trim(flux_units))
    if (CS%id_tr_adx(m) > 0) call safe_alloc_ptr(CS%tr_adx(m)%p,IsdB,IedB,jsd,jed,nz)
    if (CS%id_tr_ady(m) > 0) call safe_alloc_ptr(CS%tr_ady(m)%p,isd,ied,JsdB,JedB,nz)
    if (CS%id_tr_dfx(m) > 0) call safe_alloc_ptr(CS%tr_dfx(m)%p,IsdB,IedB,jsd,jed,nz)
    if (CS%id_tr_dfy(m) > 0) call safe_alloc_ptr(CS%tr_dfy(m)%p,isd,ied,JsdB,JedB,nz)

!    Register the tracer for horizontal advection & diffusion.
    if ((CS%id_tr_adx(m) > 0) .or. (CS%id_tr_ady(m) > 0) .or. &
        (CS%id_tr_dfx(m) > 0) .or. (CS%id_tr_dfy(m) > 0)) &
      call add_tracer_diagnostics(name, CS%tr_Reg, CS%tr_adx(m)%p, &
                                  CS%tr_ady(m)%p,CS%tr_dfx(m)%p,CS%tr_dfy(m)%p)
  enddo

end subroutine initialize_OCMIP2_CFC

subroutine init_tracer_CFC(h, tr, name, land_val, IC_val, G, CS)
  type(ocean_grid_type),                    intent(in)  :: G
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)  :: h
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out) :: tr
  character(len=*),                         intent(in)  :: name
  real,                                     intent(in)  :: land_val, IC_val
  type(OCMIP2_CFC_CS),                      pointer     :: CS

  ! This subroutine initializes a tracer array.

  logical :: OK
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (len_trim(CS%IC_file) > 0) then
    !  Read the tracer concentrations from a netcdf file.
    if (.not.file_exists(CS%IC_file, G%Domain)) &
      call MOM_error(FATAL, "initialize_OCMIP2_CFC: Unable to open "//CS%IC_file)
    if (CS%Z_IC_file) then
      OK = tracer_Z_init(tr, h, CS%IC_file, name, G)
      if (.not.OK) then
        OK = tracer_Z_init(tr, h, CS%IC_file, trim(name), G)
        if (.not.OK) call MOM_error(FATAL,"initialize_OCMIP2_CFC: "//&
                "Unable to read "//trim(name)//" from "//&
                trim(CS%IC_file)//".")
      endif
    else
      call read_data(CS%IC_file, trim(name), tr, domain=G%Domain%mpp_domain)
    endif
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      if (G%mask2dT(i,j) < 0.5) then
        tr(i,j,k) = land_val
      else
        tr(i,j,k) = IC_val
      endif
    enddo ; enddo ; enddo
  endif

end subroutine init_tracer_CFC

subroutine OCMIP2_CFC_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, CS, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),              intent(in) :: G
  type(verticalGrid_type),            intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h_old, h_new, ea, eb
  type(forcing),                      intent(in) :: fluxes
  real,                               intent(in) :: dt
  type(OCMIP2_CFC_CS),                pointer    :: CS
  real,                             optional,intent(in)  :: evap_CFL_limit
  real,                             optional,intent(in)  :: minimum_forcing_depth
!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! CFCs are relatively simple, as they are passive tracers. with only a surface
! flux as a source.

! Arguments: h_old -  Layer thickness before entrainment, in m or kg m-2.
!  (in)      h_new -  Layer thickness after entrainment, in m or kg m-2.
!  (in)      ea - an array to which the amount of fluid entrained
!                 from the layer above during this call will be
!                 added, in m or kg m-2.
!  (in)      eb - an array to which the amount of fluid entrained
!                 from the layer below during this call will be
!                 added, in m or kg m-2.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_OCMIP2_CFC.
!
! The arguments to this subroutine are redundant in that
!     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]

  real :: b1(SZI_(G))          ! b1 and c1 are variables used by the
  real :: c1(SZI_(G),SZK_(G))  ! tridiagonal solver.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    CFC11_flux, &    ! The fluxes of CFC11 and CFC12 into the ocean, in the
    CFC12_flux       ! units of CFC concentrations times meters per second.
  real, pointer, dimension(:,:,:) :: CFC11, CFC12
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_work ! Used so that h can be modified
  integer :: i, j, k, is, ie, js, je, nz, m

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return

  CFC11 => CS%CFC11 ; CFC12 => CS%CFC12

  ! These two calls unpack the fluxes from the input arrays.
  !   The -GV%Rho0 changes the sign convention of the flux and changes the units
  ! of the flux from [Conc. m s-1] to [Conc. kg m-2 s-1].
  call extract_coupler_values(fluxes%tr_fluxes, CS%ind_cfc_11_flux, ind_flux, &
                              CFC11_flux, is, ie, js, je, -GV%Rho0)
  call extract_coupler_values(fluxes%tr_fluxes, CS%ind_cfc_12_flux, ind_flux, &
                              CFC12_flux, is, ie, js, je, -GV%Rho0)

  ! Use a tridiagonal solver to determine the concentrations after the
  ! surface source is applied and diapycnal advection and diffusion occurs.
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do k=1,nz ;do j=js,je ; do i=is,ie
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo;
    call applyTracerBoundaryFluxesInOut(G, GV, CFC11, dt, fluxes, h_work, &
        evap_CFL_limit, minimum_forcing_depth)
    call tracer_vertdiff(h_work, ea, eb, dt, CFC11, G, GV, sfc_flux=CFC11_flux)

    do k=1,nz ;do j=js,je ; do i=is,ie
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo;
    call applyTracerBoundaryFluxesInOut(G, GV, CFC12, dt, fluxes, h_work, &
        evap_CFL_limit, minimum_forcing_depth)
    call tracer_vertdiff(h_work, ea, eb, dt, CFC12, G, GV, sfc_flux=CFC12_flux)
  else
    call tracer_vertdiff(h_old, ea, eb, dt, CFC11, G, GV, sfc_flux=CFC11_flux)
    call tracer_vertdiff(h_old, ea, eb, dt, CFC12, G, GV, sfc_flux=CFC12_flux)
  endif

  ! Write out any desired diagnostics.
  if (CS%mask_tracers) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      if (h_new(i,j,k) < 1.1*GV%Angstrom) then
        CS%CFC11_aux(i,j,k) = CS%CFC11_land_val
        CS%CFC12_aux(i,j,k) = CS%CFC12_land_val
      else
        CS%CFC11_aux(i,j,k) = CFC11(i,j,k)
        CS%CFC12_aux(i,j,k) = CFC12(i,j,k)
      endif
    enddo ; enddo ; enddo
    if (CS%id_CFC11>0) call post_data(CS%id_CFC11, CS%CFC11_aux, CS%diag)
    if (CS%id_CFC12>0) call post_data(CS%id_CFC12, CS%CFC12_aux, CS%diag)
  else
    if (CS%id_CFC11>0) call post_data(CS%id_CFC11, CFC11, CS%diag)
    if (CS%id_CFC12>0) call post_data(CS%id_CFC12, CFC12, CS%diag)
  endif
  do m=1,NTR
    if (CS%id_tr_adx(m)>0) &
      call post_data(CS%id_tr_adx(m),CS%tr_adx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_ady(m)>0) &
      call post_data(CS%id_tr_ady(m),CS%tr_ady(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfx(m)>0) &
      call post_data(CS%id_tr_dfx(m),CS%tr_dfx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfy(m)>0) &
      call post_data(CS%id_tr_dfy(m),CS%tr_dfy(m)%p(:,:,:),CS%diag)
  enddo

end subroutine OCMIP2_CFC_column_physics

function OCMIP2_CFC_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G
  type(verticalGrid_type),            intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h
  real, dimension(:),                 intent(out)   :: stocks
  type(OCMIP2_CFC_CS),                pointer       :: CS
  character(len=*), dimension(:),     intent(out)   :: names
  character(len=*), dimension(:),     intent(out)   :: units
  integer, optional,                  intent(in)    :: stock_index
  integer                                           :: OCMIP2_CFC_stock
! This function calculates the mass-weighted integral of all tracer stocks,
! returning the number of stocks it has calculated.  If the stock_index
! is present, only the stock corresponding to that coded index is returned.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (out)     stocks - the mass-weighted integrated amount of each tracer,
!                     in kg times concentration units.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_OCMIP2_CFC.
!  (out)     names - the names of the stocks calculated.
!  (out)     units - the units of the stocks calculated.
!  (in,opt)  stock_index - the coded index of a specific stock being sought.
! Return value: the number of stocks calculated here.

  real :: mass
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  OCMIP2_CFC_stock = 0
  if (.not.associated(CS)) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  call query_vardesc(CS%CFC11_desc, name=names(1), units=units(1), caller="OCMIP2_CFC_stock")
  call query_vardesc(CS%CFC12_desc, name=names(2), units=units(2), caller="OCMIP2_CFC_stock")
  units(1) = trim(units(1))//" kg" ; units(2) = trim(units(2))//" kg"

  stocks(1) = 0.0 ; stocks(2) = 0.0
  do k=1,nz ; do j=js,je ; do i=is,ie
    mass = G%mask2dT(i,j) * G%areaT(i,j) * h(i,j,k)
    stocks(1) = stocks(1) + CS%CFC11(i,j,k) * mass
    stocks(2) = stocks(2) + CS%CFC12(i,j,k) * mass
  enddo ; enddo ; enddo
  stocks(1) = GV%H_to_kg_m2 * stocks(1)
  stocks(2) = GV%H_to_kg_m2 * stocks(2)

  OCMIP2_CFC_stock = 2

end function OCMIP2_CFC_stock

subroutine OCMIP2_CFC_surface_state(state, h, G, CS)
  type(ocean_grid_type),                    intent(in) :: G
  type(surface),                            intent(inout) :: state
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
  type(OCMIP2_CFC_CS),                      pointer    :: CS
!   This subroutine sets up the fields that the coupler needs to calculate the
! CFC fluxes between the ocean and atmosphere.
! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_OCMIP2_CFC.

  real, dimension(SZI_(G),SZJ_(G)) :: &
    CFC11_Csurf, &  ! The CFC-11 and CFC-12 surface concentrations times the
    CFC12_Csurf, &  ! Schmidt number term, both in mol m-3.
    CFC11_alpha, &  ! The CFC-11 solubility in mol m-3 pptv-1.
    CFC12_alpha     ! The CFC-12 solubility in mol m-3 pptv-1.
  real :: ta        ! Absolute sea surface temperature in units of dekaKelvin!?!
  real :: sal       ! Surface salinity in PSU.
  real :: SST       ! Sea surface temperature in degrees Celsius.
  real :: alpha_11  ! The solubility of CFC 11 in mol m-3 pptv-1.
  real :: alpha_12  ! The solubility of CFC 12 in mol m-3 pptv-1.
  real :: sc_11, sc_12 ! The Schmidt numbers of CFC 11 and CFC 12.
  real :: sc_no_term   ! A term related to the Schmidt number.
  integer :: i, j, k, is, ie, js, je, m

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.not.associated(CS)) return

  do j=js,je ; do i=is,ie
    ta = max(0.01, (state%SST(i,j) + 273.15) * 0.01) ! Why is this in hectoKelvin?
    sal = state%SSS(i,j) ; SST = state%SST(i,j)
    !    Calculate solubilities using Warner and Weiss (1985) DSR, vol 32.
    ! The final result is in mol/cm3/pptv (1 part per trillion 1e-12)
    ! Use Bullister and Wisegavger for CCl4.
    ! The factor 1.e-09 converts from mol/(l * atm) to mol/(m3 * pptv).
    alpha_11 = exp(CS%d1_11 + CS%d2_11/ta + CS%d3_11*log(ta) + CS%d4_11*ta**2 +&
                   sal * ((CS%e3_11 * ta + CS%e2_11) * ta + CS%e1_11)) * &
               1.0e-09 * G%mask2dT(i,j)
    alpha_12 = exp(CS%d1_12 + CS%d2_12/ta + CS%d3_12*log(ta) + CS%d4_12*ta**2 +&
                   sal * ((CS%e3_12 * ta + CS%e2_12) * ta + CS%e1_12)) * &
               1.0e-09 * G%mask2dT(i,j)
    !   Calculate Schmidt numbers using coefficients given by
    ! Zheng et al (1998), JGR vol 103, C1.
    sc_11 = CS%a1_11 + SST * (CS%a2_11 + SST * (CS%a3_11 + SST * CS%a4_11)) * &
            G%mask2dT(i,j)
    sc_12 = CS%a1_12 + SST * (CS%a2_12 + SST * (CS%a3_12 + SST * CS%a4_12)) * &
            G%mask2dT(i,j)
    ! The abs here is to avoid NaNs. The model should be failing at this point.
    sc_no_term = sqrt(660.0 / (abs(sc_11) + 1.0e-30))
    CFC11_alpha(i,j) = alpha_11 * sc_no_term
    CFC11_Csurf(i,j) = CS%CFC11(i,j,1) * sc_no_term

    sc_no_term = sqrt(660.0 / (abs(sc_12) + 1.0e-30))
    CFC12_alpha(i,j) = alpha_12 * sc_no_term
    CFC12_Csurf(i,j) = CS%CFC12(i,j,1) * sc_no_term
  enddo ; enddo

  !   These calls load these values into the appropriate arrays in the
  ! coupler-type structure.
  call set_coupler_values(CFC11_alpha, state%tr_fields, CS%ind_cfc_11_flux, &
                          ind_alpha, is, ie, js, je)
  call set_coupler_values(CFC11_Csurf, state%tr_fields, CS%ind_cfc_11_flux, &
                          ind_csurf, is, ie, js, je)
  call set_coupler_values(CFC12_alpha, state%tr_fields, CS%ind_cfc_12_flux, &
                          ind_alpha, is, ie, js, je)
  call set_coupler_values(CFC12_Csurf, state%tr_fields, CS%ind_cfc_12_flux, &
                          ind_csurf, is, ie, js, je)

end subroutine OCMIP2_CFC_surface_state

subroutine OCMIP2_CFC_end(CS)
  type(OCMIP2_CFC_CS), pointer :: CS
!   This subroutine deallocates the memory owned by this module.
! Argument: CS - The control structure returned by a previous call to
!                register_OCMIP2_CFC.
  integer :: m

  if (associated(CS)) then
    if (associated(CS%CFC11)) deallocate(CS%CFC11)
    if (associated(CS%CFC12)) deallocate(CS%CFC12)
    if (associated(CS%CFC11_aux)) deallocate(CS%CFC11_aux)
    if (associated(CS%CFC12_aux)) deallocate(CS%CFC12_aux)
    do m=1,NTR
      if (associated(CS%tr_adx(m)%p)) deallocate(CS%tr_adx(m)%p)
      if (associated(CS%tr_ady(m)%p)) deallocate(CS%tr_ady(m)%p)
      if (associated(CS%tr_dfx(m)%p)) deallocate(CS%tr_dfx(m)%p)
      if (associated(CS%tr_dfy(m)%p)) deallocate(CS%tr_dfy(m)%p)
    enddo

    deallocate(CS)
  endif
end subroutine OCMIP2_CFC_end

end module MOM_OCMIP2_CFC
