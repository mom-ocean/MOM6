!> Interface to CVMix interior shear schemes
module MOM_cvmix_shear

! This file is part of MOM6. See LICENSE.md for the license.

!---------------------------------------------------
! module MOM_cvmix_shear
! Author: Brandon Reichl
! Date: Aug 31, 2016
! Purpose: Interface to CVMix interior shear schemes
! Further information to be added at a later time.
!---------------------------------------------------

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
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

public Calculate_cvmix_shear, cvmix_shear_init, cvmix_shear_is_used

!> Control structure including parameters for CVMix interior shear schemes.
type, public :: CVMix_shear_CS
  logical :: use_LMD94, use_PP81            !< Flags for various schemes
  real    :: Ri_zero                        !< LMD94 critical Richardson number
  real    :: Nu_zero                        !< LMD94 maximum interior diffusivity
  real    :: KPP_exp                        !<
  real, allocatable, dimension(:,:,:) :: N2 !< Squared Brunt-Vaisala frequency (1/s2)
  real, allocatable, dimension(:,:,:) :: S2 !< Squared shear frequency (1/s2)
  character(10) :: Mix_Scheme               !< Mixing scheme name (string)
end type CVMix_shear_CS

character(len=40)  :: mod = "MOM_CVMix_shear"  !< This module's name.

contains

!> Subroutine for calculating (internal) diffusivity
subroutine Calculate_cvmix_shear(u_H, v_H, h, tv, KH,  &
                                 KM, G, GV, CS )
  type(ocean_grid_type),                      intent(in)  :: G !< Grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV !< Vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)  :: u_H !< Initial zonal velocity on T points, in m s-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)  :: v_H !< Initial meridional velocity on T points, in m s-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)  :: h !< Layer thickness, in m or kg m-2.
  type(thermo_var_ptrs),                      intent(in)  :: tv !< Thermodynamics structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(out) :: KH !< The vertical viscosity at each interface
                                                                !! (not layer!) in m2 s-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(out) :: KM !< The vertical viscosity at each interface
                                                                !! (not layer!) in m2 s-1.
  type(CVMix_shear_CS),                       pointer     :: CS !< The control structure returned by a previous call to
                                                                !! CVMix_shear_init.
  ! Local variables
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
      Ri_Grad(:)=1.e8 !Initialize w/ large Richardson value
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
      call calculate_density(Temp_1D, Salt_1D, pres_1D, rho_1D, 1, 2*G%ke, TV%EQN_OF_STATE)

      ! N2 (can be negative) on interface
      do k = 1, G%ke
        km1 = max(1, k-1)
        kk = 2*(k-1)
        DU = (u_h(i,j,k))-(u_h(i,j,km1))
        DV = (v_h(i,j,k))-(v_h(i,j,km1))
        DRHO = (GoRho * (rho_1D(kk+1) - rho_1D(kk+2)) )
        DZ = ((0.5*(h(i,j,km1) + h(i,j,k))+GV%H_subroundoff)*GV%H_to_m)
        N2 = DRHO/DZ
        S2 = (DU*DU+DV*DV)/(DZ*DZ)
        Ri_Grad(k) = max(0.,N2)/max(S2,1.e-16)
      enddo

      ! Call to CVMix wrapper for computing interior mixing coefficients.
      call  cvmix_coeffs_shear(Mdiff_out=KM(i,j,:), &
                                   Tdiff_out=KH(i,j,:), &
                                   RICH=Ri_Grad, &
                                   nlev=G%ke,    &
                                   max_nlev=G%ke)
    enddo
  enddo

end subroutine Calculate_cvmix_shear


!> Initialized the cvmix internal shear mixing routine.
!! \note *This is where we test to make sure multiple internal shear
!!       mixing routines (including JHL) are not enabled at the same time.
!! (returns) cvmix_shear_init - True if module is to be used, False otherwise
logical function cvmix_shear_init(Time, G, GV, param_file, diag, CS)
  type(time_type),         intent(in)    :: Time !< The current time.
  type(ocean_grid_type),   intent(in)    :: G !< Grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< Vertical grid structure.
  type(param_file_type),   intent(in)    :: param_file !< Run-time parameter file handle
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure.
  type(CVMix_shear_CS),    pointer       :: CS !< This module's control structure.
  ! Local variables
  integer :: NumberTrue=0
  logical :: use_JHL
! This include declares and sets the variable "version".
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
  ! After testing for interior schemes, make sure only 0 or 1 are enabled.
  ! Otherwise, warn user and kill job.
  if ((NumberTrue).gt.1) then
     call MOM_error(FATAL, 'MOM_cvmix_shear_init: '// &
           'Multiple shear driven internal mixing schemes selected,'//&
           ' please disable all but one scheme to proceed.')
  endif
  cvmix_shear_init=(CS%use_PP81.or.CS%use_LMD94)

! Forego remainder of initialization if not using this scheme
  if (.not. cvmix_shear_init) return
  call get_param(param_file, mod, "NU_ZERO", CS%Nu_Zero, &
                 "Leading coefficient in KPP shear mixing.", &
                 units="nondim", default=5.e-3)
  call get_param(param_file, mod, "RI_ZERO", CS%Ri_Zero, &
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
  ! Allocation and initialization
  allocate( CS%N2( SZI_(G), SZJ_(G), SZK_(G)+1 ) );CS%N2(:,:,:) = 0.
  allocate( CS%S2( SZI_(G), SZJ_(G), SZK_(G)+1 ) );CS%S2(:,:,:) = 0.

end function cvmix_shear_init

!> Reads the parameters "LMD94" and "PP81" and returns state.
!!   This function allows other modules to know whether this parameterization will
!! be used without needing to duplicate the log entry.
logical function cvmix_shear_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< Run-time parameter files handle.
  ! Local variables
  logical :: LMD94, PP81
  call get_param(param_file, mod, "USE_LMD94", LMD94, &
       default=.false., do_not_log = .true.)
  call get_param(param_file, mod, "Use_PP81", PP81, &
       default=.false., do_not_log = .true.)
  cvmix_shear_is_used = (LMD94 .or. PP81)
end function cvmix_shear_is_used

end module MOM_cvmix_shear
