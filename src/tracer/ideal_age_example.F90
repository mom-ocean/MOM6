module ideal_age_example
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
!*  By Robert Hallberg, 2002                                           *
!*                                                                     *
!*    This file contains an example of the code that is needed to set  *
!*  up and use a set (in this case two) of dynamically passive tracers *
!*  for diagnostic purposes.  The tracers here are an ideal age tracer *
!*  that ages at a rate of 1/year once it is isolated from the surface,*
!*  and a vintage tracer, whose surface concentration grows exponen-   *
!*  with time with a 30-year timescale (similar to CFCs).              *
!*                                                                     *
!*    A single subroutine is called from within each file to register  *
!*  each of the tracers for reinitialization and advection and to      *
!*  register the subroutine that initializes the tracers and set up    *
!*  their output and the subroutine that does any tracer physics or    *
!*  chemistry along with diapycnal mixing (included here because some  *
!*  tracers may float or swim vertically or dye diapycnal processes).  *
!*                                                                     *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, tr                                    *
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
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
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

use coupler_util, only : set_coupler_values, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

public register_ideal_age_tracer, initialize_ideal_age_tracer
public ideal_age_tracer_column_physics, ideal_age_tracer_surface_state
public ideal_age_stock, ideal_age_example_end

! NTR_MAX is the maximum number of tracers in this module.
integer, parameter :: NTR_MAX = 3

type p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d

type, public :: ideal_age_tracer_CS ; private
  integer :: ntr    ! The number of tracers that are actually used.
  logical :: coupled_tracers = .false.  ! These tracers are not offered to the
                                        ! coupler.
  integer :: nkml       ! The number of layers in the mixed layer.  The ideal
                        ! age tracers are reset in the top nkml layers.
  character(len=200) :: IC_file ! The file in which the age-tracer initial values
                    ! can be found, or an empty string for internal initialization.
  logical :: Z_IC_file ! If true, the IC_file is in Z-space.  The default is false.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL()
  real, pointer :: tr(:,:,:,:) => NULL()   ! The array of tracers used in this
                                           ! subroutine, in g m-3?
  real, pointer :: tr_aux(:,:,:,:) => NULL() ! The masked tracer concentration
                                             ! for output, in g m-3.
  type(p3d), dimension(NTR_MAX) :: &
    tr_adx, &! Tracer zonal advective fluxes in g m-3 m3 s-1.
    tr_ady, &! Tracer meridional advective fluxes in g m-3 m3 s-1.
    tr_dfx, &! Tracer zonal diffusive fluxes in g m-3 m3 s-1.
    tr_dfy   ! Tracer meridional diffusive fluxes in g m-3 m3 s-1.
  real, dimension(NTR_MAX) :: &
    IC_val = 0.0, &    ! The (uniform) initial condition value.
    young_val = 0.0, & ! The value assigned to tr at the surface.
    land_val = -1.0, & ! The value of tr used where land is masked out.
    sfc_growth_rate, & ! The exponential growth rate for the surface value,
                       ! in units of year-1.
    tracer_start_year  ! The year in which tracers start aging, or at which the
                       ! surface value equals young_val, in years.
  logical :: mask_tracers  ! If true, tracers are masked out in massless layers.
  logical :: tracers_may_reinit  ! If true, tracers may go through the
                           ! initialization code if they are not found in the
                           ! restart files.
  logical :: tracer_ages(NTR_MAX)

  integer, dimension(NTR_MAX) :: &
    ind_tr, &  ! Indices returned by aof_set_coupler_flux if it is used and the
               ! surface tracer concentrations are to be provided to the coupler.
    id_tracer = -1, id_tr_adx = -1, id_tr_ady = -1, &
    id_tr_dfx = -1, id_tr_dfy = -1

  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()

  type(vardesc) :: tr_desc(NTR_MAX)
end type ideal_age_tracer_CS

contains

function register_ideal_age_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI
  type(verticalGrid_type),    intent(in) :: GV
  type(param_file_type),      intent(in) :: param_file
  type(ideal_age_tracer_CS),  pointer    :: CS
  type(tracer_registry_type), pointer    :: tr_Reg
  type(MOM_restart_CS),       pointer    :: restart_CS
! This subroutine is used to register tracer fields and subroutines
! to be used with MOM.
! Arguments: HI - A horizontal index type structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in/out)  tr_Reg - A pointer that is set to point to the control structure
!                  for the tracer advection and diffusion module.
!  (in)      restart_CS - A pointer to the restart control structure.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "ideal_age_example" ! This module's name.
  character(len=200) :: inputdir ! The directory where the input files are.
  character(len=48)  :: var_name ! The variable's name.
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_ideal_age_tracer
  logical :: do_ideal_age, do_vintage, do_ideal_age_dated
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "register_ideal_age_tracer called with an "// &
                             "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "DO_IDEAL_AGE", do_ideal_age, &
                 "If true, use an ideal age tracer that is set to 0 age \n"//&
                 "in the mixed layer and ages at unit rate in the interior.", &
                 default=.true.)
  call get_param(param_file, mod, "DO_IDEAL_VINTAGE", do_vintage, &
                 "If true, use an ideal vintage tracer that is set to an \n"//&
                 "exponentially increasing value in the mixed layer and \n"//&
                 "is conserved thereafter.", default=.false.)
  call get_param(param_file, mod, "DO_IDEAL_AGE_DATED", do_ideal_age_dated, &
                 "If true, use an ideal age tracer that is everywhere 0 \n"//&
                 "before IDEAL_AGE_DATED_START_YEAR, but the behaves like \n"//&
                 "the standard ideal age tracer - i.e. is set to 0 age in \n"//&
                 "the mixed layer and ages at unit rate in the interior.", &
                 default=.false.)


  call get_param(param_file, mod, "AGE_IC_FILE", CS%IC_file, &
                 "The file in which the age-tracer initial values can be \n"//&
                 "found, or an empty string for internal initialization.", &
                 default=" ")
  if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file,'/') == 0)) then
    ! Add the directory if CS%IC_file is not already a complete path.
    call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    call log_param(param_file, mod, "INPUTDIR/AGE_IC_FILE", CS%IC_file)
  endif
  call get_param(param_file, mod, "AGE_IC_FILE_IS_Z", CS%Z_IC_file, &
                 "If true, AGE_IC_FILE is in depth space, not layer space", &
                 default=.false.)
  call get_param(param_file, mod, "MASK_MASSLESS_TRACERS", CS%mask_tracers, &
                 "If true, the tracers are masked out in massless layer. \n"//&
                 "This can be a problem with time-averages.", default=.false.)
  call get_param(param_file, mod, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code \n"//&
                 "if they are not found in the restart files.  Otherwise \n"//&
                 "it is a fatal error if the tracers are not found in the \n"//&
                 "restart files of a restarted run.", default=.false.)

  CS%ntr = 0
  if (do_ideal_age) then
    CS%ntr = CS%ntr + 1 ; m = CS%ntr
    CS%tr_desc(m) = var_desc("age", "years", "Ideal Age Tracer", cmor_field_name="agessc", caller=mod)
    CS%tracer_ages(m) = .true. ; CS%sfc_growth_rate(m) = 0.0
    CS%IC_val(m) = 0.0 ; CS%young_val(m) = 0.0 ; CS%tracer_start_year(m) = 0.0
  endif

  if (do_vintage) then
    CS%ntr = CS%ntr + 1 ; m = CS%ntr
    CS%tr_desc(m) = var_desc("vintage", "years", "Exponential Vintage Tracer", &
                            caller=mod)
    CS%tracer_ages(m) = .false. ; CS%sfc_growth_rate(m) = 1.0/30.0
    CS%IC_val(m) = 0.0 ; CS%young_val(m) = 1e-20 ; CS%tracer_start_year(m) = 0.0
    call get_param(param_file, mod, "IDEAL_VINTAGE_START_YEAR", CS%tracer_start_year(m), &
                 "The date at which the ideal vintage tracer starts.", &
                 units="years", default=0.0)
  endif

  if (do_ideal_age_dated) then
    CS%ntr = CS%ntr + 1 ; m = CS%ntr
    CS%tr_desc(m) = var_desc("age_dated","years","Ideal Age Tracer with a Start Date",&
                            caller=mod)
    CS%tracer_ages(m) = .true. ; CS%sfc_growth_rate(m) = 0.0
    CS%IC_val(m) = 0.0 ; CS%young_val(m) = 0.0 ; CS%tracer_start_year(m) = 0.0
    call get_param(param_file, mod, "IDEAL_AGE_DATED_START_YEAR", CS%tracer_start_year(m), &
                 "The date at which the dated ideal age tracer starts.", &
                 units="years", default=0.0)
  endif

  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr)) ; CS%tr(:,:,:,:) = 0.0
  if (CS%mask_tracers) then
    allocate(CS%tr_aux(isd:ied,jsd:jed,nz,CS%ntr)) ; CS%tr_aux(:,:,:,:) = 0.0
  endif

  do m=1,CS%ntr
    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    call query_vardesc(CS%tr_desc(m), name=var_name, &
                       caller="register_ideal_age_tracer")
    ! Register the tracer for the restart file.
    call register_restart_field(tr_ptr, CS%tr_desc(m), &
                                .not.CS%tracers_may_reinit, restart_CS)
    ! Register the tracer for horizontal advection & diffusion.
    call register_tracer(tr_ptr, CS%tr_desc(m), param_file, HI, GV, tr_Reg, &
                         tr_desc_ptr=CS%tr_desc(m))

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = aof_set_coupler_flux(trim(var_name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_ideal_age_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_ideal_age_tracer = .true.
end function register_ideal_age_tracer

subroutine initialize_ideal_age_tracer(restart, day, G, GV, h, diag, OBC, CS, &
                                       sponge_CSp, diag_to_Z_CSp)
  logical,                            intent(in) :: restart
  type(time_type), target,            intent(in) :: day
  type(ocean_grid_type),              intent(in) :: G
  type(verticalGrid_type),            intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
  type(diag_ctrl), target,            intent(in) :: diag
  type(ocean_OBC_type),               pointer    :: OBC
  type(ideal_age_tracer_CS),          pointer    :: CS
  type(sponge_CS),                    pointer    :: sponge_CSp
  type(diag_to_Z_CS),                 pointer    :: diag_to_Z_CSp
!   This subroutine initializes the CS%ntr tracer fields in tr(:,:,:,:)
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
!                 register_ideal_age_tracer.
!  (in/out)  sponge_CSp - A pointer to the control structure for the sponges, if
!                         they are in use.  Otherwise this may be unassociated.
!  (in/out)  diag_to_Z_Csp - A pointer to the control structure for diagnostics
!                            in depth space.
  character(len=24) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for age tracer fluxes, either
                                ! years m3 s-1 or years kg s-1.
  logical :: OK
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%Time => day
  CS%diag => diag
  CS%nkml = max(GV%nkml,1)

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=name, &
                       caller="initialize_ideal_age_tracer")
    if ((.not.restart) .or. (CS%tracers_may_reinit .and. .not. &
        query_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp))) then

      if (len_trim(CS%IC_file) > 0) then
  !  Read the tracer concentrations from a netcdf file.
        if (.not.file_exists(CS%IC_file, G%Domain)) &
          call MOM_error(FATAL, "initialize_ideal_age_tracer: "// &
                                 "Unable to open "//CS%IC_file)

        if (CS%Z_IC_file) then
          OK = tracer_Z_init(CS%tr(:,:,:,m), h, CS%IC_file, name,&
                             G, -1e34, 0.0) ! CS%land_val(m))
          if (.not.OK) then
            OK = tracer_Z_init(CS%tr(:,:,:,m), h, CS%IC_file, &
                     trim(name), G, -1e34, 0.0) ! CS%land_val(m))
            if (.not.OK) call MOM_error(FATAL,"initialize_ideal_age_tracer: "//&
                    "Unable to read "//trim(name)//" from "//&
                    trim(CS%IC_file)//".")
          endif
        else
          call read_data(CS%IC_file, trim(name), CS%tr(:,:,:,m), &
                         domain=G%Domain%mpp_domain)
        endif
      else
        do k=1,nz ; do j=js,je ; do i=is,ie
          if (G%mask2dT(i,j) < 0.5) then
            CS%tr(i,j,k,m) = CS%land_val(m)
          else
            CS%tr(i,j,k,m) = CS%IC_val(m)
          endif
        enddo ; enddo ; enddo
      endif

    endif ! restart
  enddo ! Tracer loop

  if (associated(OBC)) then
  ! All tracers but the first have 0 concentration in their inflows. As this
  ! is the default value, the following calls are unnecessary.
  ! do m=1,CS%ntr
  !  call add_tracer_OBC_values(trim(CS%tr_desc(m)%name), CS%tr_Reg, 0.0)
  ! enddo
  endif

  ! This needs to be changed if the units of tracer are changed above.
  if (GV%Boussinesq) then ; flux_units = "years m3 s-1"
  else ; flux_units = "years kg s-1" ; endif

  do m=1,CS%ntr
    ! Register the tracer for the restart file.
    call query_vardesc(CS%tr_desc(m), name, units=units, longname=longname, &
                       caller="initialize_ideal_age_tracer")
    CS%id_tracer(m) = register_diag_field("ocean_model", trim(name), CS%diag%axesTL, &
        day, trim(longname) , trim(units))
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

    call register_Z_tracer(CS%tr(:,:,:,m), trim(name), longname, units, &
                           day, G, diag_to_Z_CSp)
  enddo

end subroutine initialize_ideal_age_tracer

subroutine ideal_age_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, CS, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),              intent(in) :: G
  type(verticalGrid_type),            intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h_old, h_new, ea, eb
  type(forcing),                      intent(in) :: fluxes
  real,                               intent(in) :: dt
  type(ideal_age_tracer_CS),          pointer    :: CS
  real,                             optional,intent(in)  :: evap_CFL_limit
  real,                             optional,intent(in)  :: minimum_forcing_depth
!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! This is a simple example of a set of advected passive tracers.

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
!                 register_ideal_age_tracer.
!
! The arguments to this subroutine are redundant in that
!     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_work ! Used so that h can be modified
  real :: sfc_val  ! The surface value for the tracers.
  real :: Isecs_per_year  ! The number of seconds in a year.
  real :: year            ! The time in years.
  integer :: secs, days   ! Integer components of the time type.
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,CS%ntr
      do k=1,nz ;do j=js,je ; do i=is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo;
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m) , dt, fluxes, h_work, &
          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  else
    do m=1,CS%ntr
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

  Isecs_per_year = 1.0 / (365.0*86400.0)
  !   Set the surface value of tracer 1 to increase exponentially
  ! with a 30 year time scale.
  call get_time(CS%Time, secs, days)
  year = (86400.0*days + real(secs)) * Isecs_per_year

  do m=1,CS%ntr
    if (CS%sfc_growth_rate(m) == 0.0) then
      sfc_val = CS%young_val(m)
    else
      sfc_val = CS%young_val(m) * &
          exp((year-CS%tracer_start_year(m)) * CS%sfc_growth_rate(m))
    endif
    do k=1,CS%nkml ; do j=js,je ; do i=is,ie
      if (G%mask2dT(i,j) > 0.5) then
        CS%tr(i,j,k,m) = sfc_val
      else
        CS%tr(i,j,k,m) = CS%land_val(m)
      endif
    enddo ; enddo ; enddo
  enddo
  do m=1,CS%ntr ; if (CS%tracer_ages(m) .and. &
                      (year>=CS%tracer_start_year(m))) then
!$OMP parallel do default(none) shared(is,ie,js,je,CS,nz,G,dt,Isecs_per_year,m)
    do k=CS%nkml+1,nz ; do j=js,je ; do i=is,ie
      CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + G%mask2dT(i,j)*dt*Isecs_per_year
    enddo ; enddo ; enddo
  endif ; enddo

  if (CS%mask_tracers) then
    do m=1,CS%ntr ; if (CS%id_tracer(m) > 0) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        if (h_new(i,j,k) < 1.1*GV%Angstrom) then
          CS%tr_aux(i,j,k,m) = CS%land_val(m)
        else
          CS%tr_aux(i,j,k,m) = CS%tr(i,j,k,m)
        endif
      enddo ; enddo ; enddo
    endif ; enddo
  endif

  do m=1,CS%ntr
    if (CS%mask_tracers) then
      if (CS%id_tracer(m)>0) &
        call post_data(CS%id_tracer(m),CS%tr_aux(:,:,:,m),CS%diag)
    else
      if (CS%id_tracer(m)>0) &
        call post_data(CS%id_tracer(m),CS%tr(:,:,:,m),CS%diag)
    endif
    if (CS%id_tr_adx(m)>0) &
      call post_data(CS%id_tr_adx(m),CS%tr_adx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_ady(m)>0) &
      call post_data(CS%id_tr_ady(m),CS%tr_ady(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfx(m)>0) &
      call post_data(CS%id_tr_dfx(m),CS%tr_dfx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfy(m)>0) &
      call post_data(CS%id_tr_dfy(m),CS%tr_dfy(m)%p(:,:,:),CS%diag)
  enddo

end subroutine ideal_age_tracer_column_physics

function ideal_age_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h
  real, dimension(:),                 intent(out)   :: stocks
  type(verticalGrid_type),            intent(in)    :: GV
  type(ideal_age_tracer_CS),          pointer       :: CS
  character(len=*), dimension(:),     intent(out)   :: names
  character(len=*), dimension(:),     intent(out)   :: units
  integer, optional,                  intent(in)    :: stock_index
  integer                                           :: ideal_age_stock
! This function calculates the mass-weighted integral of all tracer stocks,
! returning the number of stocks it has calculated.  If the stock_index
! is present, only the stock corresponding to that coded index is returned.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (out)     stocks - the mass-weighted integrated amount of each tracer,
!                     in kg times concentration units.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_ideal_age_tracer.
!  (out)     names - the names of the stocks calculated.
!  (out)     units - the units of the stocks calculated.
!  (in,opt)  stock_index - the coded index of a specific stock being sought.
! Return value: the number of stocks calculated here.

  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  ideal_age_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="ideal_age_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      stocks(m) = stocks(m) + CS%tr(i,j,k,m) * &
                             (G%mask2dT(i,j) * G%areaT(i,j) * h(i,j,k))
    enddo ; enddo ; enddo
    stocks(m) = GV%H_to_kg_m2 * stocks(m)
  enddo
  ideal_age_stock = CS%ntr

end function ideal_age_stock

subroutine ideal_age_tracer_surface_state(state, h, G, CS)
  type(ocean_grid_type),                    intent(in)    :: G
  type(surface),                            intent(inout) :: state
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h
  type(ideal_age_tracer_CS),                pointer       :: CS
!   This particular tracer package does not report anything back to the coupler.
! The code that is here is just a rough guide for packages that would.
! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_ideal_age_tracer.
  integer :: i, j, m, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,CS%ntr
      !   This call loads the surface vlues into the appropriate array in the
      ! coupler-type structure.
      call set_coupler_values(CS%tr(:,:,1,m), state%tr_fields, CS%ind_tr(m), &
                              ind_csurf, is, ie, js, je)
    enddo
  endif

end subroutine ideal_age_tracer_surface_state

subroutine ideal_age_example_end(CS)
  type(ideal_age_tracer_CS), pointer :: CS
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    if (associated(CS%tr_aux)) deallocate(CS%tr_aux)
    do m=1,CS%ntr
      if (associated(CS%tr_adx(m)%p)) deallocate(CS%tr_adx(m)%p)
      if (associated(CS%tr_ady(m)%p)) deallocate(CS%tr_ady(m)%p)
      if (associated(CS%tr_dfx(m)%p)) deallocate(CS%tr_dfx(m)%p)
      if (associated(CS%tr_dfy(m)%p)) deallocate(CS%tr_dfy(m)%p)
    enddo

    deallocate(CS)
  endif
end subroutine ideal_age_example_end

end module ideal_age_example
