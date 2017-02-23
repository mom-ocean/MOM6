module MOM_int_tide_input
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
!*  By Robert Hallberg, January 2013                                   *
!*                                                                     *
!*    This file contains the subroutines that sets the energy input    *
!*  to the internal tides.                                             *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, buoy, ustar, T, S, Kd, ea, eb, etc.   *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_diag_mediator, only : safe_alloc_ptr, post_data, register_diag_field
use MOM_diag_to_Z, only : diag_to_Z_CS, register_Zint_diag, calc_Zint_diags
use MOM_debugging, only : hchksum
use MOM_error_handler, only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_io, only : slasher, vardesc
use MOM_thickness_diffuse, only : vert_fill_TS
use MOM_variables, only : thermo_var_ptrs, vertvisc_type, p3d
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs

use fms_mod, only : read_data

implicit none ; private

#include <MOM_memory.h>

public set_int_tide_input, int_tide_input_init, int_tide_input_end

type, public :: int_tide_input_CS ; private
  logical :: debug           ! If true, write verbose checksums for debugging.
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  real :: TKE_itide_max ! Maximum Internal tide conversion (W m-2)
                        ! available to mix above the BBL

  real, allocatable, dimension(:,:) :: &
    TKE_itidal_coef     ! The time-invariant field that enters the TKE_itidal
                        ! input calculation, in J m-2.

  integer :: id_TKE_itidal = -1, id_Nb = -1, id_N2_bot = -1
  character(len=200) :: inputdir
end type int_tide_input_CS

type, public :: int_tide_input_type
  real, allocatable, dimension(:,:) :: &
    TKE_itidal_input, & ! The internal tide TKE input at the bottom of
                        ! the ocean, in W m-2.
    h2, &               ! The squared topographic roughness height, in m2.
    tideamp, &          ! The amplitude of the tidal velocities, in m s-1.
    Nb                  ! The bottom stratification, in s-1.
end type int_tide_input_type

contains

subroutine set_int_tide_input(u, v, h, tv, fluxes, itide, dt, G, GV, CS)
  type(ocean_grid_type),                     intent(in)    :: G
  type(verticalGrid_type),                   intent(in)    :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h
  type(thermo_var_ptrs),                     intent(in)    :: tv
  type(forcing),                             intent(in)    :: fluxes
  type(int_tide_input_type),                 intent(inout) :: itide
  real,                                      intent(in)    :: dt
  type(int_tide_input_CS),                   pointer       :: CS

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      fluxes - A structure of surface fluxes that may be used.
!  (inout)   itide - A structure containing fields related to the internal
!                    tide sources.
!  (in)      dt - The time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - This module's control structure.

  real, dimension(SZI_(G),SZJ_(G)) :: &
    N2_bot        ! The bottom squared buoyancy frequency, in s-2.

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    T_f, S_f      ! The temperature and salinity in C and PSU with the values in
                  ! the massless layers filled vertically by diffusion.
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  integer :: i, j, k, is, ie, js, je, nz
  integer :: isd, ied, jsd, jed

  real :: kappa_fill  ! diffusivity used to fill massless layers
  real :: dt_fill     ! timestep used to fill massless layers

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) call MOM_error(FATAL,"set_diffusivity: "//&
         "Module must be initialized before it is used.")

  kappa_fill = 1.e-3 ! m2 s-1
  dt_fill = 7200.

  use_EOS = associated(tv%eqn_of_state)

  ! Smooth the properties through massless layers.
  if (use_EOS) then
    call vert_fill_TS(h, tv%T, tv%S, kappa_fill, dt_fill, T_f, S_f, G, GV)
  endif

  call find_N2_bottom(h, tv, T_f, S_f, itide%h2, fluxes, G, GV, N2_bot)

!$OMP parallel do default(none) shared(is,ie,js,je,G,itide,N2_bot,CS)
  do j=js,je ; do i=is,ie
    itide%Nb(i,j) = G%mask2dT(i,j) * sqrt(N2_bot(i,j))
    itide%TKE_itidal_input(i,j) = min(CS%TKE_itidal_coef(i,j)*itide%Nb(i,j),CS%TKE_itide_max)
  enddo ; enddo

  if (CS%debug) then
    call hchksum(N2_bot,"N2_bot",G%HI,haloshift=0)
    call hchksum(itide%TKE_itidal_input,"TKE_itidal_input",G%HI,haloshift=0)
  endif

  if (CS%id_TKE_itidal > 0) call post_data(CS%id_TKE_itidal, itide%TKE_itidal_input, CS%diag)
  if (CS%id_Nb > 0) call post_data(CS%id_Nb, itide%Nb, CS%diag)
  if (CS%id_N2_bot > 0 ) call post_data(CS%id_N2_bot,N2_bot,CS%diag)

end subroutine set_int_tide_input

subroutine find_N2_bottom(h, tv, T_f, S_f, h2, fluxes, G, GV, N2_bot)
  type(ocean_grid_type),                    intent(in)   :: G
  type(verticalGrid_type),                  intent(in)   :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)   :: h
  type(thermo_var_ptrs),                    intent(in)   :: tv
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)   :: T_f, S_f
  real, dimension(SZI_(G),SZJ_(G)),         intent(in)   :: h2
  type(forcing),                            intent(in)   :: fluxes
  type(int_tide_input_CS),                  pointer      :: CS
  real, dimension(SZI_(G),SZJ_(G)),         intent(out)  :: N2_bot

  real, dimension(SZI_(G),SZK_(G)+1) :: &
    dRho_int      ! The unfiltered density differences across interfaces.
  real, dimension(SZI_(G)) :: &
    pres, &       ! The pressure at each interface, in Pa.
    Temp_int, &   ! The temperature at each interface, in degC.
    Salin_int, &  ! The salinity at each interface, in PSU.
    drho_bot, &
    h_amp, &
    hb, &
    z_from_bot, &
    dRho_dT, &    ! The partial derivatives of density with temperature and
    dRho_dS       ! salinity, in kg m-3 degC-1 and kg m-3 PSU-1.

  real :: dz_int  ! The thickness associated with an interface, in m.
  real :: G_Rho0  ! The gravitation acceleration divided by the Boussinesq
                  ! density, in m4 s-2 kg-1.
  logical :: do_i(SZI_(G)), do_any
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  G_Rho0 = GV%g_Earth / GV%Rho0

  ! Find the (limited) density jump across each interface.
  do i=is,ie
    dRho_int(i,1) = 0.0 ; dRho_int(i,nz+1) = 0.0
  enddo
!$OMP parallel do default(none) shared(is,ie,js,je,nz,tv,fluxes,G,GV,h,T_f,S_f, &
!$OMP                                  h2,N2_bot,G_Rho0) &
!$OMP                          private(pres,Temp_Int,Salin_Int,dRho_dT,dRho_dS, &
!$OMP                                  hb,dRho_bot,z_from_bot,do_i,h_amp,       &
!$OMP                                  do_any,dz_int) &
!$OMP                     firstprivate(dRho_int)
  do j=js,je
    if (associated(tv%eqn_of_state)) then
      if (associated(fluxes%p_surf)) then
        do i=is,ie ; pres(i) = fluxes%p_surf(i,j) ; enddo
      else
        do i=is,ie ; pres(i) = 0.0 ; enddo
      endif
      do K=2,nz
        do i=is,ie
          pres(i) = pres(i) + GV%H_to_Pa*h(i,j,k-1)
          Temp_Int(i) = 0.5 * (T_f(i,j,k) + T_f(i,j,k-1))
          Salin_Int(i) = 0.5 * (S_f(i,j,k) + S_f(i,j,k-1))
        enddo
        call calculate_density_derivs(Temp_int, Salin_int, pres, &
                 dRho_dT(:), dRho_dS(:), is, ie-is+1, tv%eqn_of_state)
        do i=is,ie
          dRho_int(i,K) = max(dRho_dT(i)*(T_f(i,j,k) - T_f(i,j,k-1)) + &
                              dRho_dS(i)*(S_f(i,j,k) - S_f(i,j,k-1)), 0.0)
        enddo
      enddo
    else
      do K=2,nz ; do i=is,ie
        dRho_int(i,K) = GV%Rlay(k) - GV%Rlay(k-1)
      enddo ; enddo
    endif

    ! Find the bottom boundary layer stratification.
    do i=is,ie
      hb(i) = 0.0 ; dRho_bot(i) = 0.0
      z_from_bot(i) = 0.5*GV%H_to_m*h(i,j,nz)
      do_i(i) = (G%mask2dT(i,j) > 0.5)
      h_amp(i) = sqrt(h2(i,j))
    enddo

    do k=nz,2,-1
      do_any = .false.
      do i=is,ie ; if (do_i(i)) then
        dz_int = 0.5*GV%H_to_m*(h(i,j,k) + h(i,j,k-1))
        z_from_bot(i) = z_from_bot(i) + dz_int ! middle of the layer above

        hb(i) = hb(i) + dz_int
        dRho_bot(i) = dRho_bot(i) + dRho_int(i,K)

        if (z_from_bot(i) > h_amp(i)) then
          if (k>2) then
            ! Always include at least one full layer.
            hb(i) = hb(i) + 0.5*GV%H_to_m*(h(i,j,k-1) + h(i,j,k-2))
            dRho_bot(i) = dRho_bot(i) + dRho_int(i,K-1)
          endif
          do_i(i) = .false.
        else
          do_any = .true.
        endif
      endif ; enddo
      if (.not.do_any) exit
    enddo

    do i=is,ie
      if (hb(i) > 0.0) then
        N2_bot(i,j) = (G_Rho0 * dRho_bot(i)) / hb(i)
      else ;  N2_bot(i,j) = 0.0 ; endif
    enddo
  enddo

end subroutine find_N2_bottom

subroutine int_tide_input_init(Time, G, GV, param_file, diag, CS, itide)
  type(time_type),          intent(in)    :: Time
  type(ocean_grid_type),    intent(in)    :: G
  type(verticalGrid_type),  intent(in)    :: GV
  type(param_file_type),    intent(in)    :: param_file
  type(diag_ctrl), target,  intent(inout) :: diag
  type(int_tide_input_CS),   pointer      :: CS
  type(int_tide_input_type), pointer      :: itide
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      diag_to_Z_CSp - A pointer to the Z-diagnostics control structure.
  type(vardesc) :: vd
  logical :: read_tideamp
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_int_tide_input"  ! This module's name.
  character(len=20)  :: tmpstr
  character(len=200) :: filename, tideamp_file, h2_file

  real :: mask_itidal
  real :: utide              ! constant tidal amplitude (m s-1) to be used if
                             ! tidal amplitude file is not present.
  real :: kappa_h2_factor    ! factor for the product of wavenumber * rms sgs height.
  real :: kappa_itides       ! topographic wavenumber and non-dimensional scaling
  real :: min_zbot_itides    ! Minimum ocean depth for internal tide conversion,
                             ! in m.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed

  if (associated(CS)) then
    call MOM_error(WARNING, "int_tide_input_init called with an associated "// &
                            "control structure.")
    return
  endif
  if (associated(itide)) then
    call MOM_error(WARNING, "int_tide_input_init called with an associated "// &
                            "internal tide input type.")
    return
  endif
  allocate(CS)
  allocate(itide)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")

  call get_param(param_file, mod, "INPUTDIR", CS%inputdir, default=".")
  CS%inputdir = slasher(CS%inputdir)

  call get_param(param_file, mod, "DEBUG", CS%debug, default=.false., do_not_log=.true.)

  call get_param(param_file, mod, "MIN_ZBOT_ITIDES", min_zbot_itides, &
               "Turn off internal tidal dissipation when the total \n"//&
               "ocean depth is less than this value.", units="m", default=0.0)

  call get_param(param_file, mod, "UTIDE", utide, &
               "The constant tidal amplitude used with INT_TIDE_DISSIPATION.", &
               units="m s-1", default=0.0)

  allocate(itide%Nb(isd:ied,jsd:jed))  ; itide%Nb(:,:) = 0.0
  allocate(itide%h2(isd:ied,jsd:jed))  ; itide%h2(:,:) = 0.0
  allocate(itide%TKE_itidal_input(isd:ied,jsd:jed)) ; itide%TKE_itidal_input(:,:) = 0.0
  allocate(itide%tideamp(isd:ied,jsd:jed)) ; itide%tideamp(:,:) = utide
  allocate(CS%TKE_itidal_coef(isd:ied,jsd:jed)) ; CS%TKE_itidal_coef(:,:) = 0.0

  call get_param(param_file, mod, "KAPPA_ITIDES", kappa_itides, &
               "A topographic wavenumber used with INT_TIDE_DISSIPATION. \n"//&
               "The default is 2pi/10 km, as in St.Laurent et al. 2002.", &
               units="m-1", default=8.e-4*atan(1.0))

  call get_param(param_file, mod, "KAPPA_H2_FACTOR", kappa_h2_factor, &
               "A scaling factor for the roughness amplitude with n"//&
               "INT_TIDE_DISSIPATION.",  units="nondim", default=1.0)
  call get_param(param_file, mod, "TKE_ITIDE_MAX", CS%TKE_itide_max, &
               "The maximum internal tide energy source availble to mix \n"//&
               "above the bottom boundary layer with INT_TIDE_DISSIPATION.", &
               units="W m-2",  default=1.0e3)

  call get_param(param_file, mod, "READ_TIDEAMP", read_tideamp, &
               "If true, read a file (given by TIDEAMP_FILE) containing \n"//&
               "the tidal amplitude with INT_TIDE_DISSIPATION.", default=.false.)
  if (read_tideamp) then
    call get_param(param_file, mod, "TIDEAMP_FILE", tideamp_file, &
               "The path to the file containing the spatially varying \n"//&
               "tidal amplitudes with INT_TIDE_DISSIPATION.", default="tideamp.nc")
    filename = trim(CS%inputdir) // trim(tideamp_file)
    call log_param(param_file, mod, "INPUTDIR/TIDEAMP_FILE", filename)
    call read_data(filename, 'tideamp', itide%tideamp, &
                   domain=G%domain%mpp_domain, timelevel=1)
  endif

  call get_param(param_file, mod, "H2_FILE", h2_file, &
               "The path to the file containing the sub-grid-scale \n"//&
               "topographic roughness amplitude with INT_TIDE_DISSIPATION.", &
               fail_if_missing=.true.)
  filename = trim(CS%inputdir) // trim(h2_file)
  call log_param(param_file, mod, "INPUTDIR/H2_FILE", filename)
  call read_data(filename, 'h2', itide%h2, domain=G%domain%mpp_domain, &
                 timelevel=1)

  do j=js,je ; do i=is,ie
    mask_itidal = 1.0
    if (G%bathyT(i,j) < min_zbot_itides) mask_itidal = 0.0

    itide%tideamp(i,j) = itide%tideamp(i,j) * mask_itidal * G%mask2dT(i,j)

    ! Restrict rms topo to 10 percent of column depth.
    itide%h2(i,j) = min(0.01*G%bathyT(i,j)**2, itide%h2(i,j))

    ! Compute the fixed part of internal tidal forcing; units are [J m-2] here.
    CS%TKE_itidal_coef(i,j) = 0.5*kappa_h2_factor*GV%Rho0*&
         kappa_itides * itide%h2(i,j) * itide%tideamp(i,j)**2
  enddo; enddo


  CS%id_TKE_itidal = register_diag_field('ocean_model','TKE_itidal_itide',diag%axesT1,Time, &
      'Internal Tide Driven Turbulent Kinetic Energy', 'Watt meter-2')

  CS%id_Nb = register_diag_field('ocean_model','Nb_itide',diag%axesT1,Time, &
       'Bottom Buoyancy Frequency', 'sec-1')

  CS%id_N2_bot = register_diag_field('ocean_model','N2_b_itide',diag%axesT1,Time, &
       'Bottom Buoyancy frequency squared', 's-2')

end subroutine int_tide_input_init

subroutine int_tide_input_end(CS)
  type(int_tide_input_CS), pointer :: CS

  if (associated(CS)) deallocate(CS)

end subroutine int_tide_input_end

end module MOM_int_tide_input
