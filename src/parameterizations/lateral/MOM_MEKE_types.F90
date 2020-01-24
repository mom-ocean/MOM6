module MOM_MEKE_types

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

!> This type is used to exchange information related to the MEKE calculations.
type, public :: MEKE_type
  ! Variables
  real, dimension(:,:), pointer :: &
    MEKE => NULL(), &   !< Vertically averaged eddy kinetic energy [L2 T-2 ~> m2 s-2].
    GM_src => NULL(), & !< MEKE source due to thickness mixing (GM) [R Z L2 T-3 ~> W m-2].
    mom_src => NULL(),& !< MEKE source from lateral friction in the momentum equations [R Z L2 T-3 ~> W m-2].
    GME_snk => NULL(),& !< MEKE sink from GME backscatter in the momentum equations [R Z L2 T-3 ~> W m-2].
    Kh => NULL(), &     !< The MEKE-derived lateral mixing coefficient [L2 T-1 ~> m2 s-1].
    Kh_diff => NULL(), & !< Uses the non-MEKE-derived thickness diffusion coefficient to diffuse
                        !! MEKE [L2 T-1 ~> m2 s-1].
    Rd_dx_h => NULL()   !< The deformation radius compared with the grid spacing [nondim].
                        !! Rd_dx_h is copied from VarMix_CS.
  real, dimension(:,:), pointer :: Ku => NULL() !< The MEKE-derived lateral viscosity coefficient
                        !! [L2 T-1 ~> m2 s-1]. This viscosity can be negative when representing
                        !! backscatter from unresolved eddies (see Jansen and Held, 2014).
  real, dimension(:,:), pointer :: Au => NULL() !< The MEKE-derived lateral biharmonic viscosity
                        !! coefficient [L4 T-1 ~> m4 s-1].
  ! Parameters
  real :: KhTh_fac = 1.0 !< Multiplier to map Kh(MEKE) to KhTh [nondim]
  real :: KhTr_fac = 1.0 !< Multiplier to map Kh(MEKE) to KhTr [nondim].
  real :: backscatter_Ro_pow = 0.0 !< Power in Rossby number function for backscatter.
  real :: backscatter_Ro_c = 0.0 !< Coefficient in Rossby number function for backscatter.

end type MEKE_type

end module MOM_MEKE_types
