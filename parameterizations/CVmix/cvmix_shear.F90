 module cvmix_shear

!BOP
!\newpage
! !MODULE: cvmix_shear
!
! !AUTHOR:
!  Michael N. Levy, NCAR (mlevy@ucar.edu)
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  shear mixing, and to set the viscosity and diffusivity coefficients.
!\\
!\\
!  References:\\
!  * RC Pacanowski and SGH Philander.
!  Parameterizations of Vertical Mixing in Numerical Models of Tropical Oceans.
!  Journal of Physical Oceanography, 1981.\\
!  * WG Large, JC McWilliams, and SC Doney.
!  Oceanic Vertical Mixing: A Review and a Model with a Nonlocal Boundary Layer
!  Parameterization.
!  Review of Geophysics, 1994.
!\\
!\\

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                                 &
                                    cvmix_zero,                               &
                                    cvmix_one,                                &
                                    cvmix_strlen,                             &
                                    cvmix_data_type,                          &
                                    CVMIX_OVERWRITE_OLD_VAL,                  &
                                    CVMIX_SUM_OLD_AND_NEW_VALS,               &
                                    CVMIX_MAX_OLD_AND_NEW_VALS
  use cvmix_put_get,         only : cvmix_put
  use cvmix_utils,           only : cvmix_update_wrap
!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

  public :: cvmix_init_shear
  public :: cvmix_coeffs_shear
  public :: cvmix_put_shear
  public :: cvmix_get_shear_real
  public :: cvmix_get_shear_str

  interface cvmix_coeffs_shear
    module procedure cvmix_coeffs_shear_low
    module procedure cvmix_coeffs_shear_wrap
  end interface cvmix_coeffs_shear

  interface cvmix_put_shear
    module procedure cvmix_put_shear_int
    module procedure cvmix_put_shear_real
    module procedure cvmix_put_shear_str
  end interface cvmix_put_shear

! !PUBLIC TYPES:

  ! cvmix_shear_params_type contains the necessary parameters for shear mixing
  ! (currently Pacanowski-Philander or Large et al)
  type, public :: cvmix_shear_params_type
    private
      ! Type of shear mixing to run (PP => Pacanowski-Philander, KPP => LMD94)
      character(len=cvmix_strlen) :: mix_scheme

      ! Pacanowski - Philander parameters
      ! See Eqs. (1) and (2) in 1981 paper

      ! numerator in viscosity term (O(5e-3) in PP81; default here is 0.01)
      real(cvmix_r8) :: PP_nu_zero  ! units: m^2/s

      ! coefficient of Richardson number in denominator of visc / diff terms
      ! (5 in PP81)
      real(cvmix_r8) :: PP_alpha    ! units: unitless

      ! exponent of denominator in viscosity term (2 in PP81)
      real(cvmix_r8) :: PP_exp      ! units: unitless

      ! background coefficients for visc / diff terms
      ! (1e-4 and 1e-5, respectively, in PP81; default here is 0 for both)
      real(cvmix_r8) :: PP_nu_b     ! units: m^2/s
      real(cvmix_r8) :: PP_kappa_b  ! units: m^2/s

      ! Large et al parameters
      ! See Eq. (28b) in 1994 paper

      ! leading coefficient of shear mixing formula (5e-3 in LMD94)
      real(cvmix_r8) :: KPP_nu_zero ! units: m^2/s

      ! critical Richardson number value (0.7 in LMD94)
      real(cvmix_r8) :: KPP_Ri_zero ! units: unitless

      ! Exponent of unitless factor of diffusities (3 in LMD94)
      real(cvmix_r8) :: KPP_exp     ! units: unitless

      ! Flag for what to do with old values of CVmix_vars%[MTS]diff
      integer :: handle_old_vals
  end type cvmix_shear_params_type
!EOP

  type(cvmix_shear_params_type), target :: CVmix_shear_params_saved

contains

!BOP

! !IROUTINE: cvmix_init_shear
! !INTERFACE:

  subroutine cvmix_init_shear(CVmix_shear_params_user, mix_scheme,            &
                              PP_nu_zero, PP_alpha, PP_exp, PP_nu_b,          &
                              PP_kappa_b, KPP_nu_zero, KPP_Ri_zero, KPP_exp,  &
                              old_vals)

! !DESCRIPTION:
!  Initialization routine for shear (Richardson number-based) mixing. There are
!  currently two supported schemes - set \verb|mix_scheme = 'PP'| to use the
!  Pacanowski-Philander mixing scheme or set \verb|mix_scheme = 'KPP'| to use
!  the interior mixing scheme laid out in Large et al.
!\\
!\\
!  PP requires setting $\nu_0$ (\verb|PP_nu_zero| in this routine), $alpha$
!  (\verb|PP_alpha|), and $n$ (\verb|PP_exp|), and returns
!  \begin{eqnarray*}
!  \nu_{PP} & = & \frac{\nu_0}{(1+\alpha \textrm{Ri})^n} + \nu_b \\
!  \kappa_{PP} & = & \frac{\nu}{1+\alpha \textrm{Ri}} + \kappa_b
!  \end{eqnarray*}
!  Note that $\nu_b$ and $\kappa_b$ are 0 by default, with the assumption that
!  background diffusivities are computed in the \verb|cvmix_background| module
! \\
! \\
! KPP requires setting $\nu^0$ (\verb|KPP_nu_zero|, $\textrm{Ri}_0
! ($\verb|KPP_Ri_zero|), and $p_1$ (\verb|KPP_exp|),  and returns
! $$
! \nu_{KPP} = \left\{
! \begin{array}{r l}
! \nu^0 & \textrm{Ri} < 0\\
! \nu^0 \left[1 - \frac{\textrm{Ri}}{\textrm{Ri}_0}^2\right]^{p_1}
!       & 0 < \textrm{Ri}
!           < \textrm{Ri}_0 \\
! 0     & \textrm{Ri}_0 < \textrm{Ri}
! \end{array} \right.
! $$
!
! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), optional, intent(in) :: mix_scheme,                     &
                                              old_vals
    real(cvmix_r8),   optional, intent(in) :: PP_nu_zero,                     &
                                              PP_alpha,                       &
                                              PP_exp,                         &
                                              PP_nu_b,                        &
                                              PP_kappa_b,                     &
                                              KPP_nu_zero,                    &
                                              KPP_Ri_zero,                    &
                                              KPP_exp

! !OUTPUT PARAMETERS:
    type(cvmix_shear_params_type), optional, target, intent(inout) ::         &
                                              CVmix_shear_params_user

!EOP
!BOC

    type(cvmix_shear_params_type), pointer :: CVmix_shear_params_out

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params_out => CVmix_shear_params_user
    else
      CVmix_shear_params_out => CVmix_shear_params_saved
    end if

    if (present(mix_scheme)) then
      call cvmix_put_shear("mix_scheme", trim(mix_scheme),                    &
                           CVmix_shear_params_user)
    else
      call cvmix_put_shear("mix_scheme", "KPP", CVmix_shear_params_user)
    end if

    select case (trim(CVmix_shear_params_out%mix_scheme))
      case ('PP')
        if (present(PP_nu_zero)) then
          call cvmix_put_shear("PP_nu_zero", PP_nu_zero,                      &
                               CVmix_shear_params_user)
        else
          call cvmix_put_shear("PP_nu_zero", 0.01_cvmix_r8,                   &
                               CVmix_shear_params_user)
        end if

        if (present(PP_alpha)) then
          call cvmix_put_shear("PP_alpha", PP_alpha, CVmix_shear_params_user)
        else
          call cvmix_put_shear("PP_alpha", 5, CVmix_shear_params_user)
        end if

        if (present(PP_exp)) then
          call cvmix_put_shear("PP_exp", PP_exp, CVmix_shear_params_user)
        else
          call cvmix_put_shear("PP_exp", 2, CVmix_shear_params_user)
        end if

        if (present(PP_nu_b)) then
          call cvmix_put_shear("PP_nu_b", PP_nu_b, CVmix_shear_params_user)
        else
          call cvmix_put_shear("PP_nu_b", cvmix_zero, CVmix_shear_params_user)
        end if

        if (present(PP_kappa_b)) then
          call cvmix_put_shear("PP_kappa_b", PP_kappa_b, CVmix_shear_params_user)
        else
          call cvmix_put_shear("PP_kappa_b", cvmix_zero, CVmix_shear_params_user)
        end if

      case ('KPP')
        if (present(KPP_nu_zero)) then
          call cvmix_put_shear("KPP_nu_zero", KPP_nu_zero,                    &
                               CVmix_shear_params_user)
        else
          call cvmix_put_shear("KPP_nu_zero", 50e-4_cvmix_r8,                 &
                               CVmix_shear_params_user)
        end if

        if (present(KPP_Ri_zero)) then
          call cvmix_put_shear("KPP_Ri_zero", KPP_Ri_zero,                    &
                               CVmix_shear_params_user)
        else
          call cvmix_put_shear("KPP_Ri_zero", 0.7_cvmix_r8,                   &
                               CVmix_shear_params_user)
        end if

        if (present(KPP_exp)) then
          call cvmix_put_shear("KPP_exp", KPP_exp, CVmix_shear_params_user)
        else
          call cvmix_put_shear("KPP_exp", 3, CVmix_shear_params_user)
        end if

      case DEFAULT
        print*, "ERROR: ", trim(CVmix_shear_params_out%mix_scheme),           &
                " is not a valid choice for shear mixing."
        stop 1

    end select

    if (present(old_vals)) then
      select case (trim(old_vals))
        case ("overwrite")
          call cvmix_put_shear('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,    &
                               cvmix_shear_params_user)
        case ("sum")
          call cvmix_put_shear('handle_old_vals', CVMIX_SUM_OLD_AND_NEW_VALS, &
                               cvmix_shear_params_user)
        case ("max")
          call cvmix_put_shear('handle_old_vals', CVMIX_MAX_OLD_AND_NEW_VALS, &
                               cvmix_shear_params_user)
        case DEFAULT
          print*, "ERROR: ", trim(old_vals), " is not a valid option for ",   &
                  "handling old values of diff and visc."
          stop 1
      end select
    else
      call cvmix_put_shear('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,        &
                               cvmix_shear_params_user)
    end if

!EOC

  end subroutine cvmix_init_shear

!BOP

! !IROUTINE: cvmix_coeffs_shear_wrap
! !INTERFACE:

  subroutine cvmix_coeffs_shear_wrap(CVmix_vars, CVmix_shear_params_user)

! !DESCRIPTION:
!  Computes vertical tracer and velocity mixing coefficients for
!  shear-type mixing parameterizations. Note that Richardson number
!  is needed at both T-points and U-points.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_shear_params_type), target, optional, intent(in) ::            &
                                           CVmix_shear_params_user

! !INPUT/OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars
!EOP
!BOC

    real(cvmix_r8), dimension(CVmix_vars%max_nlev+1) :: new_Mdiff, new_Tdiff
    integer :: nlev, max_nlev
    type(cvmix_shear_params_type), pointer :: CVmix_shear_params_in

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params_in => CVmix_shear_params_user
    else
      CVmix_shear_params_in => CVmix_shear_params_saved
    end if
    nlev = CVmix_vars%nlev
    max_nlev = CVmix_vars%max_nlev

    if (.not.associated(CVmix_vars%Mdiff_iface)) &
      call cvmix_put(CVmix_vars, "Mdiff", cvmix_zero, max_nlev)
    if (.not.associated(CVmix_vars%Tdiff_iface)) &
      call cvmix_put(CVmix_vars, "Tdiff", cvmix_zero, max_nlev)

    call cvmix_coeffs_shear(new_Mdiff, new_Tdiff,                             &
                            CVmix_vars%ShearRichardson_iface, nlev, max_nlev, &
                            CVmix_shear_params_user)
    call cvmix_update_wrap(CVmix_shear_params_in%handle_old_vals, max_nlev,   &
                           Mdiff_out = CVmix_vars%Mdiff_iface,                &
                           new_Mdiff = new_Mdiff,                             &
                           Tdiff_out = CVmix_vars%Tdiff_iface,                &
                           new_Tdiff = new_Tdiff)

!EOC

  end subroutine cvmix_coeffs_shear_wrap
!BOP

! !IROUTINE: cvmix_coeffs_shear_low
! !INTERFACE:

  subroutine cvmix_coeffs_shear_low(Mdiff_out, Tdiff_out, RICH, nlev,         &
                                    max_nlev, CVmix_shear_params_user)

! !DESCRIPTION:
!  Computes vertical tracer and velocity mixing coefficients for
!  shear-type mixing parameterizations. Note that Richardson number
!  is needed at both T-points and U-points.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_shear_params_type), target, optional, intent(in) ::            &
                                           CVmix_shear_params_user
    integer, intent(in) :: nlev, max_nlev
    real(cvmix_r8), dimension(max_nlev+1), intent(in) :: RICH

! !INPUT/OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(max_nlev+1), intent(inout) :: Mdiff_out,        &
                                                            Tdiff_out

!EOP
!BOC

    integer                   :: kw ! vertical cell index
    ! Parameters used in both PP81 and LMD94
    real(cvmix_r8)            :: nu_zero, loc_exp
    ! Parameters only used in PP81
    real(cvmix_r8)            :: PP_alpha, PP_nu_b, PP_kappa_b, denom
    ! Parameters only used in LMD94
    real(cvmix_r8)            :: KPP_Ri_zero
    type(cvmix_shear_params_type), pointer :: CVmix_shear_params

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params => CVmix_shear_params_user
    else
      CVmix_shear_params => CVmix_shear_params_saved
    end if

    select case (trim(CVmix_shear_params%mix_scheme))
      case ('PP')
        ! Copy parameters to make the code more legible
        nu_zero    = CVmix_shear_params%PP_nu_zero
        PP_alpha   = CVmix_shear_params%PP_alpha
        loc_exp    = CVmix_shear_params%PP_exp
        PP_nu_b    = CVmix_shear_params%PP_nu_b
        PP_kappa_b = CVmix_shear_params%PP_kappa_b

        ! Pacanowski-Philander
        do kw=1,nlev+1
          if (RICH(kw).gt.cvmix_zero) then
            denom = cvmix_one + PP_alpha * RICH(kw)
          else
            ! Treat non-negative Richardson number as Ri = 0
            denom = cvmix_one
          end if
          Mdiff_out(kw) = nu_zero / (denom**loc_exp) + PP_nu_b
          Tdiff_out(kw) = Mdiff_out(kw) / denom + PP_kappa_b
        end do

      case ('KPP')
        ! Copy parameters to make the code more legible
        nu_zero     = CVmix_shear_params%KPP_nu_zero
        KPP_Ri_zero = CVmix_shear_params%KPP_Ri_zero
        loc_exp     = CVmix_shear_params%KPP_exp

        ! Large, et al
        do kw=1,nlev+1
            if (RICH(kw).lt.cvmix_zero) then
              Tdiff_out(kw) = nu_zero
            else if (RICH(kw).lt.KPP_Ri_zero) then
              Tdiff_out(kw) = nu_zero * (cvmix_one - (RICH(kw)/KPP_Ri_zero)   &
                                                     **2)**loc_exp
            else ! Ri_g >= Ri_zero
              Tdiff_out(kw) = cvmix_zero
            end if
        end do
        ! to do: include global params for prandtl number!
        Mdiff_out = Tdiff_out

      case DEFAULT
        ! Note: this error should be caught in cvmix_init_shear
        print*, "ERROR: invalid choice for type of shear mixing."
        stop 1

    end select

!EOC

  end subroutine cvmix_coeffs_shear_low

!BOP

! !IROUTINE: cvmix_put_shear_int
! !INTERFACE:

  subroutine cvmix_put_shear_int(varname, val, CVmix_shear_params_user)

! !DESCRIPTION:
!  Write an integer value into a cvmix\_shear\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_shear_params_type), optional, target, intent(inout) ::         &
                                              CVmix_shear_params_user

!EOP
!BOC

    type(cvmix_shear_params_type), pointer :: CVmix_shear_params_out

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params_out => CVmix_shear_params_user
    else
      CVmix_shear_params_out => CVmix_shear_params_saved
    end if

    select case(trim(varname))
      case ('old_vals', 'handle_old_vals')
        CVmix_shear_params_out%handle_old_vals = val
      case DEFAULT
      call cvmix_put_shear(varname, real(val,cvmix_r8),                       &
                           CVmix_shear_params_user)
    end select

!EOC

  end subroutine cvmix_put_shear_int

!BOP

! !IROUTINE: cvmix_put_shear_real
! !INTERFACE:

  subroutine cvmix_put_shear_real(varname, val, CVmix_shear_params_user)

! !DESCRIPTION:
!  Write a real value into a cvmix\_shear\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_shear_params_type), optional, target, intent(inout) ::         &
                                              CVmix_shear_params_user

!EOP
!BOC

    type(cvmix_shear_params_type), pointer :: CVmix_shear_params_out

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params_out => CVmix_shear_params_user
    else
      CVmix_shear_params_out => CVmix_shear_params_saved
    end if

    select case (trim(varname))
      case ('PP_nu_zero')
        CVmix_shear_params_out%PP_nu_zero = val
      case ('PP_alpha')
        CVmix_shear_params_out%PP_alpha = val
      case ('PP_exp')
        CVmix_shear_params_out%PP_exp = val
      case ('PP_nu_b')
        CVmix_shear_params_out%PP_nu_b = val
      case ('PP_kappa_b')
        CVmix_shear_params_out%PP_kappa_b = val
      case ('KPP_nu_zero')
        CVmix_shear_params_out%KPP_nu_zero = val
      case ('KPP_Ri_zero')
        CVmix_shear_params_out%KPP_Ri_zero = val
      case ('KPP_exp')
        CVmix_shear_params_out%KPP_exp = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1

    end select

!EOC

  end subroutine cvmix_put_shear_real

!BOP

! !IROUTINE: cvmix_put_shear_str
! !INTERFACE:

  subroutine cvmix_put_shear_str(varname, val, CVmix_shear_params_user)

! !DESCRIPTION:
!  Write a string into a cvmix\_shear\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_shear_params_type), optional, target, intent(inout) ::         &
                                              CVmix_shear_params_user

!EOP
!BOC

    type(cvmix_shear_params_type), pointer :: CVmix_shear_params_out

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params_out => CVmix_shear_params_user
    else
      CVmix_shear_params_out => CVmix_shear_params_saved
    end if

    select case (trim(varname))
      case ('mix_scheme')
        CVmix_shear_params_out%mix_scheme = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1

    end select

!EOC

  end subroutine cvmix_put_shear_str

!BOP

! !IROUTINE: cvmix_get_shear_real
! !INTERFACE:

  function cvmix_get_shear_real(varname, CVmix_shear_params_user)

! !DESCRIPTION:
!  Read the real value of a cvmix\_shear\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),                                intent(in) :: varname
    type(cvmix_shear_params_type), optional, target, intent(in) ::            &
                                           CVmix_shear_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_get_shear_real

!EOP
!BOC

    type(cvmix_shear_params_type), pointer :: CVmix_shear_params_in

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params_in => CVmix_shear_params_user
    else
      CVmix_shear_params_in => CVmix_shear_params_saved
    end if

    cvmix_get_shear_real = cvmix_zero
    select case (trim(varname))
      case ('PP_nu_zero')
        cvmix_get_shear_real =CVmix_shear_params_in%PP_nu_zero
      case ('PP_alpha')
        cvmix_get_shear_real =CVmix_shear_params_in%PP_alpha
      case ('PP_exp')
        cvmix_get_shear_real =CVmix_shear_params_in%PP_exp
      case ('KPP_nu_zero')
        cvmix_get_shear_real =CVmix_shear_params_in%KPP_nu_zero
      case ('KPP_Ri_zero')
        cvmix_get_shear_real =CVmix_shear_params_in%KPP_Ri_zero
      case ('KPP_exp')
        cvmix_get_shear_real =CVmix_shear_params_in%KPP_exp
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1

    end select

!EOC

  end function cvmix_get_shear_real

!BOP

! !IROUTINE: cvmix_get_shear_str
! !INTERFACE:

  function cvmix_get_shear_str(varname, CVmix_shear_params_user)

! !DESCRIPTION:
!  Read the string contents of a cvmix\_shear\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),                                intent(in) :: varname
    type(cvmix_shear_params_type), optional, target, intent(in) ::            &
                                           CVmix_shear_params_user

! !OUTPUT PARAMETERS:
    character(len=cvmix_strlen) :: cvmix_get_shear_str

!EOP
!BOC

    type(cvmix_shear_params_type), pointer :: CVmix_shear_params_in

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params_in => CVmix_shear_params_user
    else
      CVmix_shear_params_in => CVmix_shear_params_saved
    end if

    select case (trim(varname))
      case ('mix_scheme')
        cvmix_get_shear_str = trim(CVmix_shear_params_in%mix_scheme)
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1

    end select

!EOC

  end function cvmix_get_shear_str


end module cvmix_shear
