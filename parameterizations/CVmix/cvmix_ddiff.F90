module cvmix_ddiff

!BOP
!\newpage
! !MODULE: cvmix_ddiff
!
! !AUTHOR:
!  Michael N. Levy, NCAR (mlevy@ucar.edu)
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  double diffusion mixing and to set the diffusivity coefficient
!  accordingly.
!\\
!\\
!  References:\\
!  * RW Schmitt.
!  Double Diffusion in Oceanography.
!  Annual Review of Fluid Mechanics, 1994.\\
!  * WG Large, JC McWilliams, and SC Doney.
!  Oceanic Vertical Mixing: A Review and a Model with a Nonlocal Boundary Layer
!  Parameterization.
!  Review of Geophysics, 1994.\\
!  * G Danabasoglu, WG Large, JJ Tribbia, PR Gent, BP Briegleb, and JC
!  McWilliams.
!  Diurnal Coupling in the Tropical Oceans of CCSM3.
!  Journal of Climate, 2006.
!\\
!\\

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                                 &
                                    cvmix_strlen,                             &
                                    cvmix_zero,                               &
                                    cvmix_one,                                &
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

  public :: cvmix_init_ddiff
  public :: cvmix_coeffs_ddiff
  public :: cvmix_put_ddiff
  public :: cvmix_get_ddiff_real

  interface cvmix_coeffs_ddiff
    module procedure cvmix_coeffs_ddiff_low
    module procedure cvmix_coeffs_ddiff_wrap
  end interface cvmix_coeffs_ddiff

  interface cvmix_put_ddiff
    module procedure cvmix_put_ddiff_str
    module procedure cvmix_put_ddiff_real
    module procedure cvmix_put_ddiff_int
  end interface cvmix_put_ddiff

! !PUBLIC TYPES:

  ! cvmix_ddiff_params_type contains the necessary parameters for double
  ! diffusion mixing
  type, public :: cvmix_ddiff_params_type
    private
      ! Max value of the stratification parameter (diffusivity = 0 for values
      ! that exceed this constant). R_p^0 in LMD94.
      real(cvmix_r8) :: strat_param_max    ! units: unitless

      ! Type of diffusive convection to use
      ! Options are Marmorino and Caldwell 1976 ("MC76"; default)
      ! and Kelley 1988, 1990 ("K90")
      character(len=cvmix_strlen) :: diff_conv_type

      ! leading coefficient in formula for salt-fingering regime for salinity
      ! diffusion (nu_f in LMD94, kappa_0 in Gokhan's paper)
      real(cvmix_r8) :: kappa_ddiff_s      ! units: m^2/s

      ! interior exponent in salt-fingering regime formula (2 in LMD94, 1 in
      ! Gokhan's paper)
      real(cvmix_r8) :: ddiff_exp1         ! units: unitless

      ! exterior exponent in salt-fingering regime formula (p2 in LMD94, 3 in
      ! Gokhan's paper)
      real(cvmix_r8) :: ddiff_exp2         ! units: unitless

      ! Exterior coefficient in diffusive convection regime (0.909 in LMD94)
      real(cvmix_r8) :: kappa_ddiff_param1 ! units: unitless

      ! Middle coefficient in diffusive convection regime (4.6 in LMD94)
      real(cvmix_r8) :: kappa_ddiff_param2 ! units: unitless

      ! Interior coefficient in diffusive convection regime (-0.54 in LMD94)
      real(cvmix_r8) :: kappa_ddiff_param3 ! units: unitless

      ! Molecular diffusivity (leading coefficient in diffusive convection
      ! regime)
      real(cvmix_r8) :: mol_diff           ! units: m^2/s

      ! Flag for what to do with old values of CVmix_vars%[MTS]diff
      integer :: handle_old_vals

  end type cvmix_ddiff_params_type

!EOP

  type(cvmix_ddiff_params_type), target :: CVmix_ddiff_params_saved

 contains

!BOP

! !IROUTINE: cvmix_init_ddiff
! !INTERFACE:

  subroutine cvmix_init_ddiff(CVmix_ddiff_params_user, strat_param_max,        &
                              kappa_ddiff_s, ddiff_exp1, ddiff_exp2, mol_diff, &
                              kappa_ddiff_param1, kappa_ddiff_param2,          &
                              kappa_ddiff_param3, diff_conv_type, old_vals)

! !DESCRIPTION:
!  Initialization routine for double diffusion mixing. This mixing technique
!  looks for two unstable cases in a column - salty water over fresher
!  water and colder water over warmer water - and computes different
!  diffusivity coefficients in each of these two locations. The parameter
!  \begin{eqnarray*}
!  R_\rho = \frac{\alpha (\partial \Theta / \partial z)}
!                {\beta (\partial S / \partial z)}
!  \end{eqnarray*}
!  to determine as a stratification parameter. If $(\partial S / \partial z)$
!  is positive and $1 < R_\rho < R_\rho^0$ then salt water sits on top
!  of fresh water and the diffusivity is given by
!  \begin{eqnarray*}
!  \kappa = \kappa^0 \left[ 1 - \left(\frac{R_\rho - 1}{R_\rho^0 - 1} \right)^{p_1}\right]^{p_2}
!  \end{eqnarray*}
!  By default, $R_\rho^0 = 2.55$, but that can be changed by setting
!  \verb|strat_param_max| in the code. Similarly, by default $p_1 = 1$
! (\verb|ddiff_exp1|), $p_2 = 3$ (\verb|ddiff_exp2|), and
!  \begin{eqnarray*}
!  \kappa^0 = \left\{ \begin{array}{r l}
!             7 \cdot 10^{-5}\ \textrm{m}^2\textrm{/s} & \textrm{for temperature}
!             \ (0.7 \cdot \verb|kappa_ddiff_s|\ \textrm{in this routine})\\
!             10^{-4}\ \textrm{m}^2\textrm{/s} & \textrm{for salinity and other tracers}
!             \ (\verb|kappa_ddiff_s|\ \textrm{in this routine}).
!                     \end{array} \right.
!  \end{eqnarray*}
!  On the other hand, if $(\partial \Theta / \partial z)$ is negative and
!  $0 < R_\rho < 1$ then cold water sits on warm warm water and the
!  diffusivity for temperature is given by
!  \begin{eqnarray*}
!  \kappa = \nu_\textrm{molecular} \cdot 0.909\exp\left\{ 4.6\exp\left[
!           -0.54\left( \frac{1}{R_\rho} - 1 \right) \right] \right\}
!  \end{eqnarray*}
!  where $\nu_\textrm{molecular}$ Is the molecular viscosity of water. By default it
!  is set to $1.5 \cdot 10^{-6}\ \textrm{m}^2\textrm{/s}$, but it can be changed
!  through \verb|mol_diff| in the code. Similarly, 0.909, 4.6, and -0.54 are the
!  default values of \verb|kappa_ddiff_param1|, \verb|kappa_ddiff_param2|, and
!  \verb|kappa_ddiff_param3|, respectively.\\
!\\
!  For salinity and other tracers, $\kappa$ above is multiplied by the factor
!  \begin{eqnarray*}
!  \textrm{factor} = \left\{ \begin{array}{c l}
!                    0.15R_\rho & R_\rho < 0.5\\
!                    1.85R_\rho - 0.85 & 0.5 \le R_\rho < 1\\
!                     \end{array} \right.
!  \end{eqnarray*}
!  $\kappa$ is stored in \verb|CVmix_vars%diff_iface(:,1)|, while the modified value
!  for non-temperature tracers is stored in \verb|CVmix_vars%diff_iface(:,2)|.
!  Note that CVMix assumes units are |'mks'|.\\
!\\
! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    real(cvmix_r8),   optional, intent(in) :: strat_param_max,                &
                                              kappa_ddiff_s,                  &
                                              ddiff_exp1,                     &
                                              ddiff_exp2,                     &
                                              mol_diff,                       &
                                              kappa_ddiff_param1,             &
                                              kappa_ddiff_param2,             &
                                              kappa_ddiff_param3
    character(len=*), optional, intent(in) :: diff_conv_type, old_vals

! !OUTPUT PARAMETERS:
    type(cvmix_ddiff_params_type), optional, target, intent(inout) ::         &
                                              CVmix_ddiff_params_user
!EOP
!BOC

    ! Unitless parameters
    if (present(strat_param_max)) then
      call cvmix_put_ddiff("strat_param_max", strat_param_max,                 &
                           CVmix_ddiff_params_user)
    else
      call cvmix_put_ddiff("strat_param_max", 2.55_cvmix_r8,                   &
                           CVmix_ddiff_params_user)
    end if

    if (present(diff_conv_type)) then
      call cvmix_put_ddiff("diff_conv_type", diff_conv_type,                  &
                            CVmix_ddiff_params_user)
    else
      call cvmix_put_ddiff("diff_conv_type", "MC76", CVmix_ddiff_params_user)
    end if

    if (present(ddiff_exp1)) then
      call cvmix_put_ddiff("ddiff_exp1", ddiff_exp1, CVmix_ddiff_params_user)
    else
      call cvmix_put_ddiff("ddiff_exp1", cvmix_one, CVmix_ddiff_params_user)
    end if

    if (present(ddiff_exp2)) then
      call cvmix_put_ddiff("ddiff_exp2", ddiff_exp2, CVmix_ddiff_params_user)
    else
      call cvmix_put_ddiff("ddiff_exp2", 3, CVmix_ddiff_params_user)
    end if

    if (present(kappa_ddiff_param1)) then
      call cvmix_put_ddiff("kappa_ddiff_param1", kappa_ddiff_param1,          &
                           CVmix_ddiff_params_user)
    else
      call cvmix_put_ddiff("kappa_ddiff_param1", 0.909_cvmix_r8,              &
                           CVmix_ddiff_params_user)
    end if

    if (present(kappa_ddiff_param2)) then
      call cvmix_put_ddiff("kappa_ddiff_param2", kappa_ddiff_param2,          &
                           CVmix_ddiff_params_user)
    else
      call cvmix_put_ddiff("kappa_ddiff_param2", 4.6_cvmix_r8,                &
                           CVmix_ddiff_params_user)
    end if

    if (present(kappa_ddiff_param3)) then
      call cvmix_put_ddiff("kappa_ddiff_param3", kappa_ddiff_param3,          &
                           CVmix_ddiff_params_user)
    else
      call cvmix_put_ddiff("kappa_ddiff_param3", -0.54_cvmix_r8,              &
                           CVmix_ddiff_params_user)
    end if

    ! Parameters with physical units
    if (present(kappa_ddiff_s)) then
      call cvmix_put_ddiff("kappa_ddiff_s", kappa_ddiff_s,                    &
                           CVmix_ddiff_params_user)
    else
      call cvmix_put_ddiff("kappa_ddiff_s", 1e-4_cvmix_r8,                    &
                           CVmix_ddiff_params_user)
    end if

    if (present(mol_diff)) then
      call cvmix_put_ddiff("mol_diff", mol_diff, CVmix_ddiff_params_user)
    else
      call cvmix_put_ddiff("mol_diff", 1.5e-6_cvmix_r8,                       &
                           CVmix_ddiff_params_user)
    end if

    if (present(old_vals)) then
      select case (trim(old_vals))
        case ("overwrite")
          call cvmix_put_ddiff('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,    &
                               cvmix_ddiff_params_user)
        case ("sum")
          call cvmix_put_ddiff('handle_old_vals', CVMIX_SUM_OLD_AND_NEW_VALS, &
                               cvmix_ddiff_params_user)
        case ("max")
          call cvmix_put_ddiff('handle_old_vals', CVMIX_MAX_OLD_AND_NEW_VALS, &
                               cvmix_ddiff_params_user)
        case DEFAULT
          print*, "ERROR: ", trim(old_vals), " is not a valid option for ",   &
                  "handling old values of diff and visc."
          stop 1
      end select
    else
      call cvmix_put_ddiff('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,        &
                               cvmix_ddiff_params_user)
    end if

!EOC

  end subroutine cvmix_init_ddiff

!BOP

! !IROUTINE: cvmix_coeffs_ddiff
! !INTERFACE:

  subroutine cvmix_coeffs_ddiff_wrap(CVmix_vars, CVmix_ddiff_params_user)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for the double diffusion mixing
!  parameterization.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_ddiff_params_type), optional, target, intent(in) ::            &
                                           CVmix_ddiff_params_user

! !INPUT/OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

    real(cvmix_r8), dimension(CVmix_vars%max_nlev+1) :: new_Tdiff, new_Sdiff
    integer :: nlev, max_nlev
    type(cvmix_ddiff_params_type), pointer :: CVmix_ddiff_params_in

    if (present(CVmix_ddiff_params_user)) then
      CVmix_ddiff_params_in => CVmix_ddiff_params_user
    else
      CVmix_ddiff_params_in => CVmix_ddiff_params_saved
    end if
    nlev = CVmix_vars%nlev
    max_nlev = CVmix_vars%max_nlev

    if (.not.associated(CVmix_vars%Tdiff_iface)) &
      call cvmix_put(CVmix_vars, "Tdiff", cvmix_zero, max_nlev)
    if (.not.associated(CVmix_vars%Sdiff_iface)) &
      call cvmix_put(CVmix_vars, "Sdiff", cvmix_zero, max_nlev)

    call cvmix_coeffs_ddiff(new_Tdiff, new_Sdiff, CVmix_vars%strat_param_num, &
                            CVmix_vars%strat_param_denom, nlev, max_nlev,     &
                            CVmix_ddiff_params_user)
    call cvmix_update_wrap(CVmix_ddiff_params_in%handle_old_vals, max_nlev,   &
                           Tdiff_out = CVmix_vars%Tdiff_iface,                &
                           new_Tdiff = new_Tdiff,                             &
                           Sdiff_out = CVmix_vars%Sdiff_iface,                &
                           new_Sdiff = new_Sdiff)

!EOC

  end subroutine cvmix_coeffs_ddiff_wrap

!BOP

! !IROUTINE: cvmix_coeffs_ddiff_low
! !INTERFACE:

  subroutine cvmix_coeffs_ddiff_low(Tdiff_out, Sdiff_out, strat_param_num,    &
                                    strat_param_denom, nlev, max_nlev,        &
                                    CVmix_ddiff_params_user)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for the double diffusion mixing
!  parameterization.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_ddiff_params_type), optional, target, intent(in) ::            &
                                           CVmix_ddiff_params_user
    integer,                             intent(in) :: nlev, max_nlev
    real(cvmix_r8), dimension(max_nlev), intent(in) :: strat_param_num,       &
                                                       strat_param_denom

! !INPUT/OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(max_nlev+1), intent(inout) :: Tdiff_out,        &
                                                            Sdiff_out

! !LOCAL VARIABLES:
    integer :: k
    real(cvmix_r8) :: ddiff, Rrho

!EOP
!BOC

    type(cvmix_ddiff_params_type), pointer :: CVmix_ddiff_params_in

    if (present(CVmix_ddiff_params_user)) then
      CVmix_ddiff_params_in => CVmix_ddiff_params_user
    else
      CVmix_ddiff_params_in => CVmix_ddiff_params_saved
    end if

    ! Determine coefficients
    Tdiff_out=cvmix_zero
    Sdiff_out=cvmix_zero
    do k = 1, nlev
      if ((strat_param_num(k).ge.strat_param_denom(k)).and.                   &
          (strat_param_denom(k).gt.cvmix_zero)) then
        ! Rrho > 1 and dS/dz < 0 => Salt fingering
        Rrho = strat_param_num(k) / strat_param_denom(k)
        if (Rrho.lt.CVmix_ddiff_params_in%strat_param_max) then
          ddiff = (cvmix_one - ((Rrho-cvmix_one)/                             &
                  (CVmix_ddiff_params_in%strat_param_max-cvmix_one))**        &
            CVmix_ddiff_params_in%ddiff_exp1)**CVmix_ddiff_params_in%ddiff_exp2
          Sdiff_out(k) = CVmix_ddiff_params_in%kappa_ddiff_s*ddiff
        end if
        Tdiff_out(k) = Sdiff_out(k)*0.7_cvmix_r8
      end if
      if ((strat_param_num(k).ge.strat_param_denom(k)).and.                   &
          (strat_param_num(k).lt.cvmix_zero)) then
        ! Rrho < 1 and dT/dz > 0 => Diffusive convection
        Rrho = strat_param_num(k) / strat_param_denom(k)
        select case (trim(CVmix_ddiff_params_in%diff_conv_type))
          case ("MC76")
            ddiff = CVmix_ddiff_params_in%mol_diff *                          &
                    CVmix_ddiff_params_in%kappa_ddiff_param1 *                &
                    exp(CVmix_ddiff_params_in%kappa_ddiff_param2*exp(         &
                        CVmix_ddiff_params_in%kappa_ddiff_param3*             &
                        (cvmix_one/Rrho-cvmix_one)))
          case ("K88")
            ddiff = CVmix_ddiff_params_in%mol_diff * 8.7_cvmix_r8 *           &
                    (Rrho**1.1_cvmix_r8)
          case DEFAULT
            print*, "ERROR: ", trim(CVmix_ddiff_params_in%diff_conv_type),    &
                    " is not a valid value for diff_conv_type"
            stop 1
        end select
        Tdiff_out(k) = ddiff
        if (Rrho.lt.0.5_cvmix_r8) then
          Sdiff_out(k) = 0.15_cvmix_r8*Rrho*ddiff
        else
          Sdiff_out(k) = (1.85_cvmix_r8*Rrho-0.85_cvmix_r8)*ddiff
        end if
      end if
    end do
    Tdiff_out(nlev+1) = cvmix_zero
    Sdiff_out(nlev+1) = cvmix_zero

!EOC

  end subroutine cvmix_coeffs_ddiff_low

!BOP

! !IROUTINE: cvmix_put_ddiff_str
! !INTERFACE:

  subroutine cvmix_put_ddiff_str(varname, val, CVmix_ddiff_params_user)

! !DESCRIPTION:
!  Write a string value into a cvmix\_ddiff\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname, val

! !OUTPUT PARAMETERS:
    type(cvmix_ddiff_params_type), optional, target, intent(inout) ::         &
                                              CVmix_ddiff_params_user
!EOP
!BOC

    type(cvmix_ddiff_params_type), pointer :: CVmix_ddiff_params_out

    if (present(CVmix_ddiff_params_user)) then
      CVmix_ddiff_params_out => CVmix_ddiff_params_user
    else
      CVmix_ddiff_params_out => CVmix_ddiff_params_saved
    end if

    select case (trim(varname))
      case ('diff_conv_type')
        select case (trim(val))
          case ('MC76')
            CVmix_ddiff_params_out%diff_conv_type = 'MC76'
          case ('K88')
            CVmix_ddiff_params_out%diff_conv_type = 'K88'
          case DEFAULT
            print*, "ERROR: ", trim(val),                                     &
                    " is not a valid value for diff_conv_type"
            stop 1
        end select

      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1

    end select

!EOC

  end subroutine cvmix_put_ddiff_str

!BOP

! !IROUTINE: cvmix_put_ddiff_real
! !INTERFACE:

  subroutine cvmix_put_ddiff_real(varname, val, CVmix_ddiff_params_user)

! !DESCRIPTION:
!  Write a real value into a cvmix\_ddiff\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_ddiff_params_type), optional, target, intent(inout) ::         &
                                              CVmix_ddiff_params_user
!EOP
!BOC

    type(cvmix_ddiff_params_type), pointer :: CVmix_ddiff_params_out

    if (present(CVmix_ddiff_params_user)) then
      CVmix_ddiff_params_out => CVmix_ddiff_params_user
    else
      CVmix_ddiff_params_out => CVmix_ddiff_params_saved
    end if

    select case (trim(varname))
      case ('strat_param_max')
        CVmix_ddiff_params_out%strat_param_max = val
      case ('ddiff_exp1')
        CVmix_ddiff_params_out%ddiff_exp1 = val
      case ('ddiff_exp2')
        CVmix_ddiff_params_out%ddiff_exp2 = val
      case ('kappa_ddiff_param1')
        CVmix_ddiff_params_out%kappa_ddiff_param1 = val
      case ('kappa_ddiff_param2')
        CVmix_ddiff_params_out%kappa_ddiff_param2 = val
      case ('kappa_ddiff_param3')
        CVmix_ddiff_params_out%kappa_ddiff_param3 = val
      case ('kappa_ddiff_s')
        CVmix_ddiff_params_out%kappa_ddiff_s = val
      case ('mol_diff')
        CVmix_ddiff_params_out%mol_diff = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1

    end select

!EOC

  end subroutine cvmix_put_ddiff_real

!BOP

! !IROUTINE: cvmix_put_ddiff_int
! !INTERFACE:

  subroutine cvmix_put_ddiff_int(varname, val, CVmix_ddiff_params_user)

! !DESCRIPTION:
!  Write an integer value into a cvmix\_ddiff\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_ddiff_params_type), optional, target, intent(inout) ::         &
                                              CVmix_ddiff_params_user

!EOP
!BOC

    type(cvmix_ddiff_params_type), pointer :: CVmix_ddiff_params_out

    if (present(CVmix_ddiff_params_user)) then
      CVmix_ddiff_params_out => CVmix_ddiff_params_user
    else
      CVmix_ddiff_params_out => CVmix_ddiff_params_saved
    end if

    select case (trim(varname))
      case ('old_vals', 'handle_old_vals')
        CVmix_ddiff_params_out%handle_old_vals = val
      case DEFAULT
        call cvmix_put_ddiff(varname, real(val,cvmix_r8),                     &
                             CVmix_ddiff_params_user)
    end select

!EOC

  end subroutine cvmix_put_ddiff_int

!BOP

! !IROUTINE: cvmix_get_ddiff_real
! !INTERFACE:

  function cvmix_get_ddiff_real(varname, CVmix_ddiff_params_user)

! !DESCRIPTION:
!  Return the real value of a cvmix\_ddiff\_params\_type variable.
!  NOTE: This function is not efficient and is only for infrequent
!  queries of ddiff parameters, such as at initialization.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),                                intent(in)    :: varname
    type(cvmix_ddiff_params_type), optional, target, intent(inout) ::         &
                                              CVmix_ddiff_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_get_ddiff_real
!EOP
!BOC

    type(cvmix_ddiff_params_type), pointer :: CVmix_ddiff_params_get

    if (present(CVmix_ddiff_params_user)) then
      CVmix_ddiff_params_get => CVmix_ddiff_params_user
    else
      CVmix_ddiff_params_get => CVmix_ddiff_params_saved
    end if

    cvmix_get_ddiff_real = cvmix_zero
    select case (trim(varname))
      case ('strat_param_max')
        cvmix_get_ddiff_real = CVmix_ddiff_params_get%strat_param_max
      case ('ddiff_exp1')
        cvmix_get_ddiff_real = CVmix_ddiff_params_get%ddiff_exp1
      case ('ddiff_exp2')
        cvmix_get_ddiff_real = CVmix_ddiff_params_get%ddiff_exp2
      case ('kappa_ddiff_param1')
        cvmix_get_ddiff_real = CVmix_ddiff_params_get%kappa_ddiff_param1
      case ('kappa_ddiff_param2')
        cvmix_get_ddiff_real = CVmix_ddiff_params_get%kappa_ddiff_param2
      case ('kappa_ddiff_param3')
        cvmix_get_ddiff_real = CVmix_ddiff_params_get%kappa_ddiff_param3
      case ('kappa_ddiff_s')
        cvmix_get_ddiff_real = CVmix_ddiff_params_get%kappa_ddiff_s
      case ('mol_diff')
        cvmix_get_ddiff_real = CVmix_ddiff_params_get%mol_diff
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1

    end select

!EOC

  end function cvmix_get_ddiff_real

end module cvmix_ddiff
