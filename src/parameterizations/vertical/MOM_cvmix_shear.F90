module MOM_cvmix_shear
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

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_checksums, only : hchksum
use MOM_error_handler, only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, EOS_type
use cvmix_shear, only : cvmix_init_shear, cvmix_coeffs_shear
use MOM_kappa_shear, only : kappa_shear_is_used
implicit none ; private

#include <MOM_memory.h>
#ifdef use_netCDF
#include <netcdf.inc>
#endif

public Calculate_cvmix_shear, cvmix_shear_init

type, public :: CVMix_shear_CS ! ; private
  logical :: use_LMD94, use_PP81
  real    :: Ri_zero
  real    :: Nu_zero
  real    :: KPP_exp
  real, allocatable, dimension(:,:,:) :: N2        !< Squared Brunt-Vaisala frequency (1/s2)
  real, allocatable, dimension(:,:,:) :: S2        !< Squared shear frequency (1/s2)
  character(10) :: Mix_Scheme
end type CVMix_shear_CS

! integer :: id_clock_project, id_clock_KQ, id_clock_avg, id_clock_setup
  character(len=40)  :: mod = "MOM_CVMix_shear"  ! This module's name.

#undef  DEBUG
#undef  ADD_DIAGNOSTICS

contains

subroutine Calculate_cvmix_shear(u_in, v_in, h, tv, KH,  &
                                 KM, G, GV, CS )
  type(ocean_grid_type),                      intent(in)    :: G
  type(verticalGrid_type),                    intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: u_in
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: v_in
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: h
  type(thermo_var_ptrs),                      intent(in)    :: tv
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(out)   :: KH
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(out)   :: KM
  type(CVMix_shear_CS),                       pointer       :: CS
!
! ----------------------------------------------
! Subroutine for calculating diffusivity (and TKE?)
! ----------------------------------------------
! Arguments: u_in - Initial zonal velocity, in m s-1. (Intent in)
!  (in)      v_in - Initial meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in/out)  KH - The diapycnal diffusivity at each interface
!                       (not layer!) in m2 s-1.  Initially this is the value
!                       from the previous timestep, which may accelerate the
!                       iteration toward convergence.
!  (in/out)  KM - The vertical viscosity at each interface
!                    (not layer!) in m2 s-1. This discards any previous value
!                    i.e. intent(out) and simply sets Kv = Prandtl * Kd_turb
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 CVMix_shear_init.

  integer :: i, j, k, kk, km1
  real :: gorho
  real :: pref, DU, DV, DRHO, DZ, N2, S2
  real, dimension(2*(G%ke)) :: pres_1d, temp_1d, salt_1d, rho_1d 
  real, dimension(G%ke+1) ::  Ri_Grad !< Gradient Richardson number

  ! some constants
  GoRho = GV%g_Earth / GV%Rho0

  do j = G%jsc, G%jec
    do i = G%isc, G%iec

      ! skip calling for land points
      if (G%mask2dT(i,j)==0.) cycle

      ! Richardson number computed for each cell in a column.
      pRef = 0.
      do k=1,G%ke

        ! pressure, temp, and saln for EOS
        ! kk+1 = k fields
        ! kk+2 = km1 fields
        km1  = max(1, k-1)
        kk   = 2*(k-1)
        pres_1D(kk+1) = pRef
        pres_1D(kk+2) = pRef
        Temp_1D(kk+1) = TV%T(i,j,k)
        Temp_1D(kk+2) = TV%T(i,j,km1)
        Salt_1D(kk+1) = TV%S(i,j,k)
        Salt_1D(kk+2) = TV%S(i,j,km1)

        ! pRef is pressure at interface between k and km1.
        ! iterate pRef for next pass through k-loop.
        pRef = pRef + GV%H_to_Pa * h(i,j,k)

      enddo ! k-loop finishes

      ! compute in-situ density
      call calculate_density(Temp_1D, Salt_1D, pres_1D, rho_1D, 1, 3*G%ke, TV%EQN_OF_STATE)

      ! N2 (can be negative) and N (non-negative) on interfaces.
      ! deltaRho is non-local rho difference used for bulk Richardson number.
      ! N_1d is local N (with floor) used for unresolved shear calculation.
      do k = 1, G%ke
        km1 = max(1, k-1)
        kk = 2*(k-1)
        DU = (u_in(i,j,k)+u_in(i-1,j,k))-(u_in(i,j,km1)+u_in(i-1,j,km1))
        DV = (v_in(i,j,k)+v_in(i,j-1,k))-(v_in(i,j,km1)+v_in(i,j-1,km1))
        DRHO = (GoRho * (rho_1D(kk+1) - rho_1D(kk+2)) )
        DZ = ((0.5*(h(i,j,km1) + h(i,j,k))+GV%H_subroundoff)*GV%H_to_m)
        N2 = DRHO/DZ
        S2 = (DU*DU+DV*DV)/(DZ*DZ) 
        Ri_Grad(G%ke) = N2/max(S2,1.e-8)
      enddo
      Ri_Grad(G%ke+1) = 1.e8

      call  cvmix_coeffs_shear(Mdiff_out=KM(i,j,:), &
                                   Tdiff_out=KH(i,j,:), & 
                                   RICH=Ri_Grad, &
                                   nlev=G%ke,    &
                                   max_nlev=G%ke)


    ENDDO;
  ENDDO;

  return

end subroutine Calculate_cvmix_shear


logical function cvmix_shear_init(Time, G, GV, param_file, diag, CS)
  !
  ! Initialized the cvmix internal shear mixing routine.
  !  *tests to make sure multiple internal shear mixing routines
  !   are not enabled at the same time.
  !
  type(time_type),         intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(verticalGrid_type), intent(in)    :: GV
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ctrl), target, intent(inout) :: diag
  type(CVMix_shear_CS),    pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (returns) cvmix_shear_init - True if module is to be used, False otherwise
! This include declares and sets the variable "version".
  INTEGER :: NumberTrue
  LOGICAL :: use_JHL
#include "version_variable.h"

  if (associated(CS)) then
    call MOM_error(WARNING, "cvmix_shear_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

! Set default, read and log parameters
  call log_version(param_file, mod, version, &
    "Parameterization of shear-driven turbulence via CVMix (various options)")
  call get_param(param_file, mod, "USE_LMD94", CS%use_LMD94, &
                 "If true, use the Large-McWilliams-Doney (JGR 1994) \n"//&
                 "shear mixing parameterization.", default=.false.)
  if (CS%use_LMD94) then
     NumberTrue=NumberTrue + 1
     CS%Mix_Scheme='KPP'
  endif
  call get_param(param_file, mod, "USE_PP81", CS%use_PP81, &
                 "If true, use the Pacanowski and Philander (JPO 1981) \n"//&
                 "shear mixing parameterization.", default=.false.)
  if (CS%use_PP81) then
     NumberTrue = NumberTrue + 1
     CS%Mix_Scheme='PP'
  endif
  use_JHL=kappa_shear_is_used(param_file)
  if (use_JHL) NumberTrue = NumberTrue + 1
  if ((NumberTrue).gt.1) then
     call MOM_error(FATAL, 'MOM_cvmix_shear_init: '// &
           'Multiple shear driven internal mixing schemes selected,'//&
           ' please disable all but one scheme to proceed.')
  endif
  cvmix_shear_init=(CS%use_PP81.or.CS%use_LMD94)
  print*,CS%use_LMD94
  stop
! Forego remainder of initialization if not using this scheme
  if (.not. cvmix_shear_init) return
  call get_param(param_file, mod, "NU_ZERO", CS%Nu_Zero, &
                 "Leading coefficient in KPP shear mixing.", &
                 units="nondim", default=5.e-3)
  call get_param(param_file, mod, "RI_ZERO", CS%Nu_Zero, &
                 "Critical Richardson for KPP shear mixing,"// &
                 " NOTE this the internal mixing and this is"// &
                 " not for setting the boundary layer depth." &
                 ,units="nondim", default=0.7)
  call get_param(param_file, mod, "KPP_EXP", CS%KPP_exp, &
                 "Exponent of unitless factor of diffusivities,"// &
                 " for KPP internal shear mixing scheme." &
                 ,units="nondim", default=3.0)
  call cvmix_init_shear(mix_scheme=CS%mix_scheme, &
                        KPP_nu_zero=CS%Nu_Zero,   &
                        KPP_Ri_zero=CS%Ri_zero,   &
                        KPP_exp=CS%KPP_exp)
  !Allocation and initialization
  allocate( CS%N2( SZI_(G), SZJ_(G), SZK_(G)+1 ) );CS%N2(:,:,:) = 0.
  allocate( CS%S2( SZI_(G), SZJ_(G), SZK_(G)+1 ) );CS%S2(:,:,:) = 0.

end function cvmix_shear_init

end module MOM_cvmix_shear
