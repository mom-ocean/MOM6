module cvmix_background

!BOP
!\newpage
! !MODULE: cvmix_background
!
! !AUTHOR:
!  Michael N. Levy, NCAR (mlevy@ucar.edu)
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  time independent static background mixing coefficients.  It specifies
!  either a scalar, 1D, or 2D field for viscosity and diffusivity. It also
!  calculates the background diffusivity using the Bryan-Lewis method.
!  It then sets the viscosity and diffusivity to the specified value.
!\\
!\\
!  References:\\
!  * K Bryan and LJ Lewis.
!  A Water Mass Model of the World Ocean.
!  Journal of Geophysical Research, 1979.
!\\
!\\

! !USES:

  use cvmix_kinds_and_types, only : cvmix_PI,                                 &
                                    cvmix_r8,                                 &
                                    cvmix_strlen,                             &
                                    cvmix_zero,                               &
                                    cvmix_data_type,                          &
                                    cvmix_global_params_type,                 &
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

  public :: cvmix_init_bkgnd
  public :: cvmix_coeffs_bkgnd
  public :: cvmix_bkgnd_lvary_horizontal
  public :: cvmix_bkgnd_static_Mdiff
  public :: cvmix_bkgnd_static_Tdiff
  public :: cvmix_put_bkgnd
  public :: cvmix_get_bkgnd_real_2D

  interface cvmix_init_bkgnd
    module procedure cvmix_init_bkgnd_scalar
    module procedure cvmix_init_bkgnd_1D
    module procedure cvmix_init_bkgnd_2D
    module procedure cvmix_init_bkgnd_BryanLewis_wrap
    module procedure cvmix_init_bkgnd_BryanLewis_low
  end interface cvmix_init_bkgnd

  interface cvmix_coeffs_bkgnd
    module procedure cvmix_coeffs_bkgnd_low
    module procedure cvmix_coeffs_bkgnd_wrap
  end interface cvmix_coeffs_bkgnd

  interface cvmix_put_bkgnd
    module procedure cvmix_put_bkgnd_int
    module procedure cvmix_put_bkgnd_real
    module procedure cvmix_put_bkgnd_real_1D
    module procedure cvmix_put_bkgnd_real_2D
  end interface cvmix_put_bkgnd

! !PUBLIC TYPES:

  ! cvmix_bkgnd_params_type contains the necessary parameters for background
  ! mixing. Background mixing fields can vary from level to level as well as
  ! over latitude and longitude.
  type, public :: cvmix_bkgnd_params_type
    private
      ! 3D viscosity field (horizontal dimensions are collapsed into first
      ! dimension, vertical is second dimension)
      real(cvmix_r8), allocatable :: static_Mdiff(:,:) ! ncol, max_nlev+1
                                                       ! units: m^2/s
      ! 3D diffusivity field (horizontal dimensions are collapsed into first
      ! dimension, vertical is second dimension)
      real(cvmix_r8), allocatable :: static_Tdiff(:,:) ! ncol, max_nlev+1
                                                       ! units: m^2/s

      ! Flag for what to do with old values of CVmix_vars%[MTS]diff
      integer :: handle_old_vals

      ! Note: need to include some logic to avoid excessive memory use
      !       when static_[MT]diff are constant or 1-D
      logical :: lvary_vertical   ! True => multiple levels
      logical :: lvary_horizontal ! True => multiple columns
  end type cvmix_bkgnd_params_type

!EOP

  type(cvmix_bkgnd_params_type), target :: CVmix_bkgnd_params_saved

contains

!BOP

! !IROUTINE: cvmix_init_bkgnd_scalar
! !INTERFACE:

  subroutine cvmix_init_bkgnd_scalar(bkgnd_Tdiff, bkgnd_Mdiff, old_vals,      &
                                     CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Initialization routine for static background mixing coefficients. For each
!  column, this routine sets the static viscosity / diffusivity to the given
!  scalar constants.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    real(cvmix_r8),             intent(in) :: bkgnd_Tdiff
    real(cvmix_r8),             intent(in) :: bkgnd_Mdiff
    character(len=*), optional, intent(in) :: old_vals

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), optional, target, intent(inout) :: &
                                              CVmix_bkgnd_params_user

!EOP
!BOC

    type(cvmix_bkgnd_params_type), pointer :: CVmix_bkgnd_params_out

    CVmix_bkgnd_params_out => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_out => CVmix_bkgnd_params_user
    end if

    ! Clean up memory in bkgnd_params_type (will be re-allocated in put call)
    if (allocated(CVmix_bkgnd_params_out%static_Mdiff))                       &
      deallocate(CVmix_bkgnd_params_out%static_Mdiff)
    if (allocated(CVmix_bkgnd_params_out%static_Tdiff))                       &
      deallocate(CVmix_bkgnd_params_out%static_Tdiff)

    ! Set static_Mdiff and static_Tdiff in background_input_type
    call cvmix_put_bkgnd('static_Mdiff', bkgnd_Mdiff, CVmix_bkgnd_params_user)
    call cvmix_put_bkgnd('static_Tdiff', bkgnd_Tdiff, CVmix_bkgnd_params_user)

    if (present(old_vals)) then
      select case (trim(old_vals))
        case ("overwrite")
          call cvmix_put_bkgnd('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,    &
                               cvmix_bkgnd_params_user)
        case ("sum")
          call cvmix_put_bkgnd('handle_old_vals', CVMIX_SUM_OLD_AND_NEW_VALS, &
                               cvmix_bkgnd_params_user)
        case ("max")
          call cvmix_put_bkgnd('handle_old_vals', CVMIX_MAX_OLD_AND_NEW_VALS, &
                               cvmix_bkgnd_params_user)
        case DEFAULT
          print*, "ERROR: ", trim(old_vals), " is not a valid option for ",   &
                  "handling old values of diff and visc."
          stop 1
      end select
    else
      call cvmix_put_bkgnd('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,        &
                               cvmix_bkgnd_params_user)
    end if

!EOC

  end subroutine cvmix_init_bkgnd_scalar

!BOP

! !IROUTINE: cvmix_init_bkgnd_1D
! !INTERFACE:

  subroutine cvmix_init_bkgnd_1D(bkgnd_Tdiff, bkgnd_Mdiff, ncol, old_vals,    &
                                 CVmix_params_user, CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Initialization routine for static background mixing coefficients. For each
!  column, this routine sets the static viscosity / diffusivity to the given
!  1D field. If field varies horizontally, need to include ncol!
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    real(cvmix_r8), dimension(:),          intent(in) :: bkgnd_Tdiff
    real(cvmix_r8), dimension(:),          intent(in) :: bkgnd_Mdiff
    integer,                     optional, intent(in) :: ncol
    character(len=cvmix_strlen), optional, intent(in) :: old_vals
    type(cvmix_global_params_type), optional, target, intent(in) :: &
                                                  CVmix_params_user

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type),  optional, target, intent(inout) :: &
                                               CVmix_bkgnd_params_user
!EOP
!BOC

    ! local vars
    integer :: nlev
    type(cvmix_global_params_type), pointer :: CVmix_params_in
    type(cvmix_bkgnd_params_type),  pointer :: CVmix_bkgnd_params_out

    nullify(CVmix_params_in)
    if (present(CVmix_params_user)) then
      CVmix_params_in => CVmix_params_user
      nlev = CVmix_params_in%max_nlev
    else
      if (.not.present(ncol)) then
        print*, "ERROR: You must specify either ncol or a global param type", &
                "containing max_nlev!"
        stop 1
      end if
    endif

    CVmix_bkgnd_params_out => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_out => CVmix_bkgnd_params_user
    end if

    ! NOTE: need to verify that bkgnd_[MT]diff are ncol x 1 or 1 x nlev+1

    ! Clean up memory in bkgnd_params_type (will be re-allocated in put call)
    if (allocated(CVmix_bkgnd_params_out%static_Mdiff))                       &
      deallocate(CVmix_bkgnd_params_out%static_Mdiff)
    if (allocated(CVmix_bkgnd_params_out%static_Tdiff))                       &
      deallocate(CVmix_bkgnd_params_out%static_Tdiff)

    ! Set static_[MT]diff in background_input_type
    if (present(ncol)) then
      call cvmix_put_bkgnd('static_Mdiff', bkgnd_Mdiff,                       &
                           CVmix_bkgnd_params_user, ncol=ncol)
      call cvmix_put_bkgnd('static_Tdiff', bkgnd_Tdiff,                       &
                           CVmix_bkgnd_params_user, ncol=ncol)
    else
      call cvmix_put_bkgnd('static_Mdiff', bkgnd_Mdiff,                       &
                           CVmix_bkgnd_params_user, nlev=nlev)
      call cvmix_put_bkgnd('static_Tdiff', bkgnd_Tdiff,                       &
                           CVmix_bkgnd_params_user, nlev=nlev)
    end if

    if (present(old_vals)) then
      select case (trim(old_vals))
        case ("overwrite")
          call cvmix_put_bkgnd('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,    &
                               cvmix_bkgnd_params_user)
        case ("sum")
          call cvmix_put_bkgnd('handle_old_vals', CVMIX_SUM_OLD_AND_NEW_VALS, &
                               cvmix_bkgnd_params_user)
        case ("max")
          call cvmix_put_bkgnd('handle_old_vals', CVMIX_MAX_OLD_AND_NEW_VALS, &
                               cvmix_bkgnd_params_user)
        case DEFAULT
          print*, "ERROR: ", trim(old_vals), " is not a valid option for ",   &
                  "handling old values of diff and visc."
          stop 1
      end select
    else
      call cvmix_put_bkgnd('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,        &
                               cvmix_bkgnd_params_user)
    end if

!EOC

  end subroutine cvmix_init_bkgnd_1D

!BOP

! !IROUTINE: cvmix_init_bkgnd_2D
! !INTERFACE:

  subroutine cvmix_init_bkgnd_2D(bkgnd_Tdiff, bkgnd_Mdiff, ncol,              &
                                 CVmix_params_in, old_vals,                   &
                                 CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Initialization routine for static background mixing coefficients. For each
!  column, this routine sets the static viscosity / diffusivity to the given
!  2D field.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    real(cvmix_r8), dimension(:,:),        intent(in) :: bkgnd_Tdiff
    real(cvmix_r8), dimension(:,:),        intent(in) :: bkgnd_Mdiff
    integer,                               intent(in) :: ncol
    character(len=cvmix_strlen), optional, intent(in) :: old_vals
    type(cvmix_global_params_type),        intent(in) :: CVmix_params_in

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type),  target, optional, intent(inout) ::        &
                                               CVmix_bkgnd_params_user
!EOP
!BOC

    ! local vars
    integer :: nlev
    type(cvmix_bkgnd_params_type),  pointer :: CVmix_bkgnd_params_out

    CVmix_bkgnd_params_out => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_out => CVmix_bkgnd_params_user
    end if

    ! NOTE: need to verify that bkgnd_[MT]diff are ncol x nlev+1

    nlev = CVmix_params_in%max_nlev

    ! Clean up memory in bkgnd_params_type (will be re-allocated in put call)
    if (allocated(CVmix_bkgnd_params_out%static_Mdiff))                       &
      deallocate(CVmix_bkgnd_params_out%static_Mdiff)
    if (allocated(CVmix_bkgnd_params_out%static_Tdiff))                       &
      deallocate(CVmix_bkgnd_params_out%static_Tdiff)

    ! Set static_[MT]diff in background_input_type
    call cvmix_put_bkgnd("static_Mdiff", bkgnd_Mdiff, ncol, nlev,             &
                         CVmix_bkgnd_params_user)
    call cvmix_put_bkgnd("static_Tdiff", bkgnd_Tdiff, ncol, nlev,             &
                         CVmix_bkgnd_params_user)

    if (present(old_vals)) then
      select case (trim(old_vals))
        case ("overwrite")
          call cvmix_put_bkgnd('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,    &
                               cvmix_bkgnd_params_user)
        case ("sum")
          call cvmix_put_bkgnd('handle_old_vals', CVMIX_SUM_OLD_AND_NEW_VALS, &
                               cvmix_bkgnd_params_user)
        case ("max")
          call cvmix_put_bkgnd('handle_old_vals', CVMIX_MAX_OLD_AND_NEW_VALS, &
                               cvmix_bkgnd_params_user)
        case DEFAULT
          print*, "ERROR: ", trim(old_vals), " is not a valid option for ",   &
                  "handling old values of diff and visc."
          stop 1
      end select
    else
      call cvmix_put_bkgnd('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,        &
                               cvmix_bkgnd_params_user)
    end if

!EOC

  end subroutine cvmix_init_bkgnd_2D

!BOP

! !IROUTINE: cvmix_init_bkgnd_BryanLewis_wrap
! !INTERFACE:

  subroutine cvmix_init_bkgnd_BryanLewis_wrap(CVmix_vars, bl1, bl2, bl3, bl4, &
                                              CVmix_params_in, old_vals,      &
                                              CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Calls cvmix_init_bkgnd_BryanLewis_low
!
! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    ! Contains depth and nlev
    type(cvmix_data_type), intent(in) :: CVmix_vars
    ! Units are first column if CVmix_data%depth is m, second if cm
    real(cvmix_r8), intent(in) :: bl1,     &! m^2/s or cm^2/s
                                  bl2,     &! m^2/s or cm^2/s
                                  bl3,     &! 1/m   or 1/cm
                                  bl4       ! m     or cm
    character(len=cvmix_strlen),            optional, intent(in) :: old_vals
    type(cvmix_global_params_type), intent(in) :: CVmix_params_in

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type),  target, optional, intent(inout) ::        &
                                               CVmix_bkgnd_params_user
!EOP
!BOC

    call cvmix_init_bkgnd(CVmix_params_in%max_nlev,-CVMix_vars%zw_iface, bl1,      &
                          bl2, bl3, bl4, CVmix_params_in%prandtl,                  &
                          old_vals, CVmix_bkgnd_params_user)

!EOC

  end subroutine cvmix_init_bkgnd_BryanLewis_wrap

!BOP

! !IROUTINE: cvmix_coeffs_bkgnd_low
! !INTERFACE:

  subroutine cvmix_init_bkgnd_BryanLewis_low(max_nlev, zw, bl1, bl2, bl3, bl4, &
                                      prandtl, old_vals, CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Initialization routine for Bryan-Lewis diffusivity/viscosity calculation.
!  For each column, this routine sets the static viscosity \& diffusivity
!  based on the specified parameters. Note that the units of these parameters
!  must be consistent with the units of viscosity and diffusivity -- either
!  cgs or mks, but do not mix and match!
!  \\
!  \\
!  The Bryan-Lewis parameterization is based on the following:
!  \begin{eqnarray*}
!  \kappa_{BL} &=& \textrm{bl1} + \frac{\textrm{bl2}}{\pi}\tan^{-1}\bigg(
!                  \textrm{bl3}(|z|-\textrm{bl4})\bigg)\\
!  \nu_{BL} &=& \textrm{Pr}\cdot\kappa_{BL}
!  \end{eqnarray*}
!  This method is based on the following paper:
!  \begin{quote}
!  \emph{A Water Mass Model of the World Ocean}\\
!  K. Bryan and L. J. Lewis\\
!  Journal of Geophysical Research, vol 84 (1979), pages 2503-2517.
!  \end{quote}
!
!  In that paper, they recommend the parameters
!  \begin{itemize}
!  \item[] bl1 $= 8 \cdot 10^{-5}$ m$^2/$s
!  \item[] bl2 $= 1.05 \cdot 10^{-4}$ m$^2/$s
!  \item[] bl3 $= 4.5 \cdot 10^{-3}$ m$^{-1}$
!  \item[] bl4 $= 2500$ m
!  \end{itemize}
!  However, more recent usage of their scheme may warrant different settings.
!
! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    integer,                      intent(in) :: max_nlev
    real(cvmix_r8), dimension(max_nlev+1), intent(in) :: zw
    ! Units are first column if CVmix_data%depth is m, second if cm
    real(cvmix_r8), intent(in) :: bl1,     &! m^2/s or cm^2/s
                                  bl2,     &! m^2/s or cm^2/s
                                  bl3,     &! 1/m   or 1/cm
                                  bl4,     &! m     or cm
                                  prandtl   ! nondim
    character(len=cvmix_strlen),            optional, intent(in) :: old_vals

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type),  target, optional, intent(inout) ::        &
                                               CVmix_bkgnd_params_user
!EOP
!BOC

    ! Pointers to parameter data type
    type(cvmix_bkgnd_params_type),  pointer :: CVmix_bkgnd_params_out

    ! Local copies to make code easier to read
    real(cvmix_r8), dimension(max_nlev+1) :: Mdiff, Tdiff

    CVmix_bkgnd_params_out => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_out => CVmix_bkgnd_params_user
    end if

    ! Clean up memory in bkgnd_params_type (will be re-allocated in put call)
    if (allocated(CVmix_bkgnd_params_out%static_Mdiff))                       &
      deallocate(CVmix_bkgnd_params_out%static_Mdiff)
    if (allocated(CVmix_bkgnd_params_out%static_Tdiff))                       &
      deallocate(CVmix_bkgnd_params_out%static_Tdiff)

    ! Set static_[MT]diff in background_input_type
    Tdiff = bl1 + (bl2/cvmix_PI)*atan(bl3*(zw-bl4))
    Mdiff = prandtl*Tdiff

    call cvmix_put_bkgnd("static_Mdiff", Mdiff, CVmix_bkgnd_params_user,      &
                         nlev=max_nlev)
    call cvmix_put_bkgnd("static_Tdiff", Tdiff, CVmix_bkgnd_params_user,      &
                         nlev=max_nlev)

    if (present(old_vals)) then
      select case (trim(old_vals))
        case ("overwrite")
          call cvmix_put_bkgnd('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,    &
                               cvmix_bkgnd_params_user)
        case ("sum")
          call cvmix_put_bkgnd('handle_old_vals', CVMIX_SUM_OLD_AND_NEW_VALS, &
                               cvmix_bkgnd_params_user)
        case ("max")
          call cvmix_put_bkgnd('handle_old_vals', CVMIX_MAX_OLD_AND_NEW_VALS, &
                               cvmix_bkgnd_params_user)
        case DEFAULT
          print*, "ERROR: ", trim(old_vals), " is not a valid option for ",   &
                  "handling old values of diff and visc."
          stop 1
      end select
    else
      call cvmix_put_bkgnd('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,        &
                               cvmix_bkgnd_params_user)
    end if

!EOC

  end subroutine cvmix_init_bkgnd_BryanLewis_low

!BOP

! !IROUTINE: cvmix_coeffs_bkgnd_wrap
! !INTERFACE:

  subroutine cvmix_coeffs_bkgnd_wrap(CVmix_vars, colid,                       &
                                     CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Computes vertical tracer and velocity mixing coefficients for static
!  background mixing. This routine simply copies viscosity / diffusivity
!  values from CVmix\_bkgnd\_params to CVmix\_vars.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:

    ! Need to know column for pulling data from static_[MT]diff
    integer,                               optional, intent(in) :: colid
    type(cvmix_bkgnd_params_type), target, optional, intent(in) ::            &
                                           CVmix_bkgnd_params_user

! !INPUT/OUTPUT PARAMETERS:

    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

    real(cvmix_r8), dimension(CVmix_vars%max_nlev+1) :: new_Mdiff, new_Tdiff
    type(cvmix_bkgnd_params_type), pointer :: CVmix_bkgnd_params_in
    integer :: nlev, max_nlev

    CVmix_bkgnd_params_in => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_in => CVmix_bkgnd_params_user
    end if

    nlev = CVmix_vars%nlev
    max_nlev = CVmix_vars%max_nlev

    if (.not.associated(CVmix_vars%Mdiff_iface)) &
      call cvmix_put(CVmix_vars, "Mdiff", cvmix_zero, max_nlev)
    if (.not.associated(CVmix_vars%Tdiff_iface)) &
      call cvmix_put(CVmix_vars, "Tdiff", cvmix_zero, max_nlev)

    call cvmix_coeffs_bkgnd(new_Mdiff, new_Tdiff, nlev, max_nlev, colid,      &
                            CVmix_bkgnd_params_user)
    call cvmix_update_wrap(CVmix_bkgnd_params_in%handle_old_vals, max_nlev,   &
                           Mdiff_out = CVmix_vars%Mdiff_iface,                &
                           new_Mdiff = new_Mdiff,                             &
                           Tdiff_out = CVmix_vars%Tdiff_iface,                &
                           new_Tdiff = new_Tdiff)

!EOC

  end subroutine cvmix_coeffs_bkgnd_wrap

!BOP

! !IROUTINE: cvmix_coeffs_bkgnd_low
! !INTERFACE:

  subroutine cvmix_coeffs_bkgnd_low(Mdiff_out, Tdiff_out, nlev, max_nlev,     &
                                    colid, CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Computes vertical tracer and velocity mixing coefficients for static
!  background mixing. This routine simply copies viscosity / diffusivity
!  values from CVmix\_bkgnd\_params to CVmix\_vars.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:

    ! Need to know column for pulling data from static_[MT]diff
    integer,                                         intent(in) :: nlev,      &
                                                                   max_nlev
    integer,                               optional, intent(in) :: colid
    type(cvmix_bkgnd_params_type), target, optional, intent(in) ::            &
                                           CVmix_bkgnd_params_user

! !OUTPUT PARAMETERS:
    ! Using intent(inout) because memory should already be allocated
    real(cvmix_r8), dimension(max_nlev+1), intent(inout) :: Mdiff_out,        &
                                                            Tdiff_out

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer :: kw

    do kw=1,nlev+1
      Mdiff_out(kw) = cvmix_bkgnd_static_Mdiff(CVmix_bkgnd_params_user, kw,   &
                                               colid)
      Tdiff_out(kw) = cvmix_bkgnd_static_Tdiff(CVmix_bkgnd_params_user, kw,   &
                                               colid)
    end do

!EOC

  end subroutine cvmix_coeffs_bkgnd_low

!BOP

! !IROUTINE: cvmix_bkgnd_lvary_horizontal
! !INTERFACE:

  function cvmix_bkgnd_lvary_horizontal(CVmix_bkgnd_params_test)

! !DESCRIPTION:
!  Returns whether the background viscosity and diffusivity are
!  varying with horizontal position.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), intent(in) :: CVmix_bkgnd_params_test

! !OUTPUT PARAMETERS:
  logical :: cvmix_bkgnd_lvary_horizontal
!EOP
!BOC

    cvmix_bkgnd_lvary_horizontal = CVmix_bkgnd_params_test%lvary_horizontal

!EOC

  end function cvmix_bkgnd_lvary_horizontal

!BOP

! !IROUTINE: cvmix_bkgnd_static_Mdiff
! !INTERFACE:

  function cvmix_bkgnd_static_Mdiff(CVmix_bkgnd_params_user,kw,colid)

! !DESCRIPTION:
!  Obtain the background diffusivity value at a position in a water column.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), target, optional, intent(in) ::            &
                                           CVmix_bkgnd_params_user
    integer, optional, intent(in) :: kw, colid

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_bkgnd_static_Mdiff
!EOP
!BOC

    type(cvmix_bkgnd_params_type),  pointer :: CVmix_bkgnd_params_in
    integer :: cid, kid

    ! Error check
    CVmix_bkgnd_params_in => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_in => CVmix_bkgnd_params_user
    end if

    if (CVmix_bkgnd_params_in%lvary_horizontal) then
      if (present(colid)) then
        cid = colid
      else
        print*, "ERROR: need to pass colid when static_Mdiff varies across",  &
                " columns."
        stop 1
      end if
    else
      cid = 1
    end if

    if (CVmix_bkgnd_params_in%lvary_vertical) then
      if (present(kw)) then
        kid = kw
      else
        print*, "ERROR: need to pass kw (level id) when static_Mdiff varies", &
                "across levels columns."
        stop 1
      end if
    else
      kid = 1
    end if

    cvmix_bkgnd_static_Mdiff = CVmix_bkgnd_params_in%static_Mdiff(cid, kid)

!EOC

  end function cvmix_bkgnd_static_Mdiff

!BOP

! !IROUTINE: cvmix_bkgnd_static_Tdiff
! !INTERFACE:

  function cvmix_bkgnd_static_Tdiff(CVmix_bkgnd_params_user,kw,colid)

! !DESCRIPTION:
!  Obtain the background diffusivity value at a position in a water column.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), target, optional, intent(in) ::            &
                                           CVmix_bkgnd_params_user
    integer, optional, intent(in) :: kw, colid

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_bkgnd_static_Tdiff
!EOP
!BOC

    type(cvmix_bkgnd_params_type),  pointer :: CVmix_bkgnd_params_in
    integer :: cid, kid

    ! Error che
    CVmix_bkgnd_params_in => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_in => CVmix_bkgnd_params_user
    end if

    if (CVmix_bkgnd_params_in%lvary_horizontal) then
      if (present(colid)) then
        cid = colid
      else
        print*, "ERROR: need to pass colid when static_Tdiff varies across",  &
                " columns."
        stop 1
      end if
    else
      cid = 1
    end if

    if (CVmix_bkgnd_params_in%lvary_vertical) then
      if (present(kw)) then
        kid = kw
      else
        print*, "ERROR: need to pass kw (level id) when static_Tdiff varies",  &
                "across levels columns."
        stop 1
      end if
    else
      kid = 1
    end if

    cvmix_bkgnd_static_Tdiff = CVmix_bkgnd_params_in%static_Tdiff(cid, kid)

!EOC

  end function cvmix_bkgnd_static_Tdiff

!BOP

! !IROUTINE: cvmix_put_bkgnd_int
! !INTERFACE:

  subroutine cvmix_put_bkgnd_int(varname, val, CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Write a real value into a cvmix\_bkgnd\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), target, optional, intent(inout) ::         &
                                              CVmix_bkgnd_params_user

!EOP
!BOC

    type(cvmix_bkgnd_params_type), pointer :: CVmix_bkgnd_params_out

    CVmix_bkgnd_params_out => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_out => CVmix_bkgnd_params_user
    end if

    select case (trim(varname))
      case ('old_vals', 'handle_old_vals')
        CVmix_bkgnd_params_out%handle_old_vals = val
      case DEFAULT
        call cvmix_put_bkgnd(varname, real(val,cvmix_r8),                     &
                             CVmix_bkgnd_params_user)
    end select

!EOC

  end subroutine cvmix_put_bkgnd_int

!BOP

! !IROUTINE: cvmix_put_bkgnd_real
! !INTERFACE:

  subroutine cvmix_put_bkgnd_real(varname, val, CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Write a real value into a cvmix\_bkgnd\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), target, optional, intent(inout) ::         &
                                              CVmix_bkgnd_params_user

!EOP
!BOC

    type(cvmix_bkgnd_params_type), pointer :: CVmix_bkgnd_params_out

    CVmix_bkgnd_params_out => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_out => CVmix_bkgnd_params_user
    end if

    select case (trim(varname))
      case ('static_Mdiff')
        if (.not.allocated(CVmix_bkgnd_params_out%static_Mdiff)) then
          allocate(CVmix_bkgnd_params_out%static_Mdiff(1,1))
          CVmix_bkgnd_params_out%lvary_horizontal=.false.
          CVmix_bkgnd_params_out%lvary_vertical=.false.
        else
          print*, "WARNING: overwriting static_Mdiff!"
        end if
        CVmix_bkgnd_params_out%static_Mdiff(:,:) = val

      case ('static_Tdiff')
        if (.not.allocated(CVmix_bkgnd_params_out%static_Tdiff)) then
          allocate(CVmix_bkgnd_params_out%static_Tdiff(1,1))
          CVmix_bkgnd_params_out%lvary_horizontal=.false.
          CVmix_bkgnd_params_out%lvary_vertical=.false.
        else
          print*, "WARNING: overwriting static_Tdiff!"
        end if
        CVmix_bkgnd_params_out%static_Tdiff(:,:) = val

      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1

    end select

!EOC

  end subroutine cvmix_put_bkgnd_real

!BOP

! !IROUTINE: cvmix_put_bkgnd_real_1D
! !INTERFACE:

  subroutine cvmix_put_bkgnd_real_1D(varname, val, CVmix_bkgnd_params_user,   &
                                    ncol, nlev)

! !DESCRIPTION:
!  Write an array of real values into a cvmix\_bkgnd\_params\_type variable.
!  You must use \verb|opt='horiz'| to specify that the field varies in the
!  horizontal direction, otherwise it is assumed to vary in the vertical.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),             intent(in) :: varname
    real(cvmix_r8), dimension(:), intent(in) :: val
    integer, optional,            intent(in) :: ncol, nlev

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), target, optional, intent(inout) ::         &
                                              CVmix_bkgnd_params_user

!EOP
!BOC

    ! Local vars
    integer, dimension(2) :: dims
    integer               :: data_dims
    logical               :: lvary_horizontal

    type(cvmix_bkgnd_params_type), pointer :: CVmix_bkgnd_params_out

    ! Error checking to make sure dimension is specified
    if ((.not.present(ncol)).and.(.not.present(nlev))) then
      print*, "ERROR: when putting 1D data in cvmix_bkgnd_params_type ", &
              "you must specify nlev or ncol!"
      stop 1
    end if

    if ((present(ncol)).and.(present(nlev))) then
      print*, "ERROR: when putting 1D data in cvmix_bkgnd_params_type ", &
              "you can not specify both nlev or ncol!"
      stop 1
    end if

    CVmix_bkgnd_params_out => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_out => CVmix_bkgnd_params_user
    end if

    data_dims = size(val)
    if (present(ncol)) then
      if (data_dims.gt.ncol) then
        print*, "ERROR: data array is bigger than number of columns specified."
        stop 1
      end if
      lvary_horizontal=.true.
      dims(1) = ncol
      dims(2) = 1
    else
      if (data_dims.gt.nlev+1) then
        print*, "ERROR: data array is bigger than number of levels specified."
        stop 1
      end if
      lvary_horizontal=.false.
      dims(1) = 1
      dims(2) = nlev+1
    end if

    select case (trim(varname))
      case ('static_Mdiff')
        if (.not.allocated(CVmix_bkgnd_params_out%static_Mdiff)) then
          allocate(CVmix_bkgnd_params_out%static_Mdiff(dims(1),dims(2)))
          CVmix_bkgnd_params_out%lvary_horizontal = lvary_horizontal
          CVmix_bkgnd_params_out%lvary_vertical = .not.lvary_horizontal
        else
          print*, "WARNING: overwriting static_Mdiff!"
        end if
        if (any(shape(CVmix_bkgnd_params_out%static_Mdiff).ne.dims)) then
          print*, "ERROR: dimensions of static_Mdiff do not match what was ", &
                  "sent to cvmix_put"
          stop 1
        end if
        if (lvary_horizontal) then
          CVmix_bkgnd_params_out%static_Mdiff(:,1)           = cvmix_zero
          CVmix_bkgnd_params_out%static_Mdiff(1:data_dims,1) = val
        else
          CVmix_bkgnd_params_out%static_Mdiff(1,:)           = cvmix_zero
          CVmix_bkgnd_params_out%static_Mdiff(1,1:data_dims) = val
        end if

      case ('static_Tdiff')
        if (.not.allocated(CVmix_bkgnd_params_out%static_Tdiff)) then
          allocate(CVmix_bkgnd_params_out%static_Tdiff(dims(1),dims(2)))
          CVmix_bkgnd_params_out%lvary_horizontal = lvary_horizontal
          CVmix_bkgnd_params_out%lvary_vertical = .not.lvary_horizontal
        else
          print*, "WARNING: overwriting static_Tdiff!"
        end if
        if (any(shape(CVmix_bkgnd_params_out%static_Tdiff).ne.dims)) then
          print*, "ERROR: dimensions of static_Tdiff do not match what was ", &
                  "sent to cvmix_put"
          stop 1
        end if
        if (lvary_horizontal) then
          CVmix_bkgnd_params_out%static_Tdiff(:,1)           = cvmix_zero
          CVmix_bkgnd_params_out%static_Tdiff(1:data_dims,1) = val
        else
          CVmix_bkgnd_params_out%static_Tdiff(1,:)           = cvmix_zero
          CVmix_bkgnd_params_out%static_Tdiff(1,1:data_dims) = val
        end if

      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1

    end select

!EOC

  end subroutine cvmix_put_bkgnd_real_1D

!BOP

! !IROUTINE: cvmix_put_bkgnd_real_2D
! !INTERFACE:

  subroutine cvmix_put_bkgnd_real_2D(varname, val, ncol, nlev,                &
                                     CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Write a 2D array of real values into a cvmix\_bkgnd\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),               intent(in) :: varname
    real(cvmix_r8), dimension(:,:), intent(in) :: val
    integer,                        intent(in) :: ncol, nlev

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type),  optional, target, intent(inout) ::        &
                                               CVmix_bkgnd_params_user

!EOP
!BOC

    ! Local vars
    integer, dimension(2) :: dims, data_dims
    type(cvmix_bkgnd_params_type),  pointer :: CVmix_bkgnd_params_out

    CVmix_bkgnd_params_out => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_out => CVmix_bkgnd_params_user
    end if

    dims      = (/ncol, nlev+1/)
    data_dims = shape(val)

    if (any(data_dims.gt.dims)) then
      print*, "ERROR: data being put in cvmix_bkgnd_params_type is larger ", &
              "than (ncol, nlev+1)"
      stop 1
    end if

    select case (trim(varname))
      case ('static_Mdiff')
        if (.not.allocated(CVmix_bkgnd_params_out%static_Mdiff)) then
          allocate(CVmix_bkgnd_params_out%static_Mdiff(dims(1),dims(2)))
          CVmix_bkgnd_params_out%lvary_horizontal=.true.
          CVmix_bkgnd_params_out%lvary_vertical=.true.
        else
          print*, "WARNING: overwriting static_Mdiff!"
        end if
        if (any(shape(CVmix_bkgnd_params_out%static_Mdiff).ne.dims)) then
          print*, "ERROR: dimensions of static_Mdiff do not match what was ", &
                  "sent to cvmix_put"
          stop 1
        end if
        CVmix_bkgnd_params_out%static_Mdiff = cvmix_zero
        CVmix_bkgnd_params_out%static_Mdiff(1:data_dims(1),1:data_dims(2))= val

      case ('static_Tdiff')
        if (.not.allocated(CVmix_bkgnd_params_out%static_Tdiff)) then
          allocate(CVmix_bkgnd_params_out%static_Tdiff(dims(1),dims(2)))
          CVmix_bkgnd_params_out%lvary_horizontal=.true.
          CVmix_bkgnd_params_out%lvary_vertical=.true.
        else
          print*, "WARNING: overwriting static_Tdiff!"
        end if
        if (any(shape(CVmix_bkgnd_params_out%static_Tdiff).ne.dims)) then
          print*, "ERROR: dimensions of static_Tdiff do not match what was ", &
                  "sent to cvmix_put"
          stop 1
        end if
        CVmix_bkgnd_params_out%static_Tdiff = cvmix_zero
        CVmix_bkgnd_params_out%static_Tdiff(1:data_dims(1),1:data_dims(2))= val

      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1

    end select

!EOC

  end subroutine cvmix_put_bkgnd_real_2D

!BOP

! !IROUTINE: cvmix_get_bkgnd_real_2D
! !INTERFACE:

  function cvmix_get_bkgnd_real_2D(varname, CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Read the real values of a cvmix\_bkgnd\_params\_type 2D array variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),                                intent(in) :: varname
    type(cvmix_bkgnd_params_type), target, optional, intent(in) ::            &
                                             CVmix_bkgnd_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8), allocatable, dimension(:,:) :: cvmix_get_bkgnd_real_2D

!EOP
!BOC

    type(cvmix_bkgnd_params_type), pointer :: CVmix_bkgnd_params_get
    integer :: dim1, dim2

    CVmix_bkgnd_params_get => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_get => CVmix_bkgnd_params_user
    end if
    dim1 = size(CVmix_bkgnd_params_get%static_Mdiff,1)
    dim2 = size(CVmix_bkgnd_params_get%static_Mdiff,2)
    allocate(cvmix_get_bkgnd_real_2D(dim1, dim2))

    select case (trim(varname))
      case ('static_Mdiff')
        cvmix_get_bkgnd_real_2D = CVmix_bkgnd_params_get%static_Mdiff(:,:)
      case ('static_Tdiff')
        cvmix_get_bkgnd_real_2D = CVmix_bkgnd_params_get%static_Tdiff(:,:)
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
    end select

!EOC

  end function cvmix_get_bkgnd_real_2D


end module cvmix_background

