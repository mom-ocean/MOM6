!> A column-wise toolbox for implementing neutral diffusion
module MOM_neutral_diffusion_aux

use MOM_EOS,                   only : EOS_type, extract_member_EOS, EOS_LINEAR, EOS_TEOS10, EOS_WRIGHT
use MOM_EOS,                   only : calculate_density_derivs, calculate_density_second_derivs
use MOM_error_handler,         only : MOM_error, FATAL, WARNING
use polynomial_functions,      only : evaluation_polynomial, first_derivative_polynomial

! This file is part of MOM6. See LICENSE.md for the license.
implicit none ; private

public mark_unstable_cells
public mark_unstable_cells_i
public calc_delta_rho
public refine_nondim_position
public check_neutral_positions
public kahan_sum

contains

!> Given the reconsturcitons of dRdT, dRdS, T, S mark the cells which are stably stratified parts of the water column
!! For an layer to be unstable the top interface must be denser than the bottom or the bottom interface of the layer
!!
subroutine mark_unstable_cells(nk, dRdT, dRdS,T, S, stable_cell, ns)
  integer,                intent(in)    :: nk          !< Number of levels in a column
  real, dimension(nk),    intent(in)    :: dRdT        !< drho/dT (kg/m3/degC)
  real, dimension(nk),    intent(in)    :: dRdS        !< drho/dS (kg/m3/ppt)
  real, dimension(nk),    intent(in)    :: T           !< drho/dS (kg/m3/ppt)
  real, dimension(nk),    intent(in)    :: S           !< drho/dS (kg/m3/ppt)
  logical, dimension(nk), intent(  out) :: stable_cell !< True if this cell is unstably stratified
  integer,                intent(  out) :: ns          !< Number of neutral surfaces in unmasked part of the column

  integer :: k, first_stable, prev_stable
  real :: delta_rho

  ns = 0
  ! If only one cell, then we really shouldn't do anything
  if (nk==1) then
    stable_cell(nk)=.true.
    ns = 2
    return
  endif
  first_stable = 1
  prev_stable = 1
  ! First sweep down and find the first place where the column is stable
  do k=1,nk-1
    delta_rho = ( (dRdT(k) + dRdT(k+1))*(T(k)-T(k+1)) ) + ( (dRdS(k) + dRdS(k+1))*(S(k)-S(k+1)) )
    if (delta_rho <= 0.) then
      first_stable = k+1
      prev_stable = k
      stable_cell(k) = .true.
      ns = ns + 2
      exit
    else
      stable_cell(k) = .false.
    endif
  enddo

  ! Loop through the rest of the column
  do k=first_stable,nk
    delta_rho = ( (dRdT(prev_stable) + dRdT(k))*(T(prev_stable)-T(k)) ) + ( (dRdS(prev_stable) + dRdS(k))*(S(prev_stable)-S(k)) )
    if (delta_rho <= 0.) then
      stable_cell(k) = .true.
      prev_stable = k
      ns = ns + 2
    else
      stable_cell(k) = .false.
    endif
  enddo

end subroutine mark_unstable_cells

subroutine mark_unstable_cells_i(nk, dRdT, dRdS,T, S, stable_cell, ns)
  integer,                intent(in)    :: nk          !< Number of levels in a column
  real, dimension(nk,2),    intent(in)    :: dRdT        !< drho/dT (kg/m3/degC)
  real, dimension(nk,2),    intent(in)    :: dRdS        !< drho/dS (kg/m3/ppt)
  real, dimension(nk,2),    intent(in)    :: T           !< drho/dS (kg/m3/ppt)
  real, dimension(nk,2),    intent(in)    :: S           !< drho/dS (kg/m3/ppt)
  logical, dimension(nk), intent(  out) :: stable_cell !< True if this cell is unstably stratified
  integer,                intent(  out) :: ns          !< Number of neutral surfaces in unmasked part of the column

  integer :: k, first_stable, prev_stable
  real :: delta_rho

  ! If only one cell, then we really shouldn't do anything
  if (nk==1) then
    stable_cell(nk)=.true.
    ns = 2
    return
  endif

  do k=1,nk
    ! Only check cell which are stable
    if (stable_cell(k)) then
      delta_rho = ( (dRdT(k,1) + dRdT(k,2))*(T(k,1)-T(k,2)) ) + ( (dRdS(k,1) + dRdS(k,2))*(S(k,1)-S(k,2)) )
      if (delta_rho > 0.) then
        stable_cell(k) = .false.
        ns = ns - 2
      endif
    endif
  enddo

end subroutine mark_unstable_cells_i

!> Calculate the difference in neutral density between a reference T, S, alpha, and beta
!! and a point on the polynomial reconstructions of T, S
subroutine calc_delta_rho(deg, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S, x0, ref_pres, EOS, &
                          delta_rho, P_out, T_out, S_out, alpha_avg_out, beta_avg_out, delta_T_out, delta_S_out)
  integer,                intent(in)  :: deg       !< Degree of polynomial reconstruction
  real,                   intent(in)  :: T_ref     !< Temperature at reference surface
  real,                   intent(in)  :: S_ref     !< Salinity at reference surface
  real,                   intent(in)  :: alpha_ref !< dRho/dT at reference surface
  real,                   intent(in)  :: beta_ref  !< dRho/dS at reference surface
  real,                   intent(in)  :: P_top     !< Pressure (Pa) at top interface of layer to be searched
  real,                   intent(in)  :: P_bot     !< Pressure (Pa) at bottom interface
  real, dimension(deg+1), intent(in)  :: ppoly_T   !< Coefficients of T reconstruction
  real, dimension(deg+1), intent(in)  :: ppoly_S   !< Coefficients of S reconstruciton
  real,                   intent(in)  :: x0        !< Nondimensional position to evaluate
  real,                   intent(in)  :: ref_pres  !< Reference pressure
  type(EOS_type),         pointer     :: EOS       !< Equation of state structure
  real,                   intent(out) :: delta_rho
  real,         optional, intent(out) :: P_out         !< Pressure at point x0
  real,         optional, intent(out) :: T_out         !< Temperature at point x0
  real,         optional, intent(out) :: S_out         !< Salinity at point x0
  real,         optional, intent(out) :: alpha_avg_out !< Average of alpha between reference and x0
  real,         optional, intent(out) :: beta_avg_out  !< Average of beta between reference and x0
  real,         optional, intent(out) :: delta_T_out   !< Difference in temperature between reference and x0
  real,         optional, intent(out) :: delta_S_out   !< Difference in salinity between reference and x0

  real :: alpha, beta, alpha_avg, beta_avg, P_int, T, S, delta_T, delta_S

  P_int = (1. - x0)*P_top + x0*P_bot
  T = evaluation_polynomial( ppoly_T, deg+1, x0 )
  S = evaluation_polynomial( ppoly_S, deg+1, x0 )
  ! Interpolated pressure if using locally referenced neutral density
  if (ref_pres<0.) then
    call calculate_density_derivs( T, S, P_int, alpha, beta, EOS )
  else
  ! Constant reference pressure (isopycnal)
    call calculate_density_derivs( T, S, ref_pres, alpha, beta, EOS )
  endif

  ! Calculate the f(P) term for Newton's method
  alpha_avg = 0.5*( alpha + alpha_ref )
  beta_avg = 0.5*( beta + beta_ref )
  delta_T = T - T_ref
  delta_S = S - S_ref
  delta_rho = alpha_avg*delta_T + beta_avg*delta_S

  ! If doing a Newton step, these quantities are needed, otherwise they can just be optional
  if (present(P_out)) P_out = P_int
  if (present(T_out)) T_out = T
  if (present(S_out)) S_out = S
  if (present(alpha_avg_out)) alpha_avg_out = alpha_avg
  if (present(beta_avg_out))  beta_avg_out = beta_avg
  if (present(delta_T_out)) delta_T_out = delta_T
  if (present(delta_S_out)) delta_S_out = delta_S

end subroutine calc_delta_rho

!> Use root-finding methods to find where dRho = 0, based on the equation of state and the polynomial
!! reconstructions of temperature, salinity. Initial guess is based on the zero crossing of based on linear
!! profiles of dRho, T, and S, between the top and bottom interface. If second derivatives of the EOS are available,
!! it starts with a Newton's method. However, Newton's method is not guaranteed to be bracketed, a check is performed
!! to see if it it diverges outside the interval. In that case (or in the case that second derivatives are not
!! available), Brent's method is used following the implementation found at
!! https://people.sc.fsu.edu/~jburkardt/f_src/brent/brent.f90
real function refine_nondim_position(max_iter, tolerance, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, deg, &
      ppoly_T, ppoly_S, EOS, x0, drho_top, drho_bot, min_bound, ref_pres, force_brent)
  integer,            intent(in) :: max_iter    !< Number of maximum iterations to use
  real,               intent(in) :: tolerance   !< Convergence criterion for delta_rho
  real,               intent(in) :: T_ref       !< Temperature of the neutral surface at the searched from interface
  real,               intent(in) :: S_ref       !< Salinity of the neutral surface at the searched from interface
  real,               intent(in) :: alpha_ref   !< dRho/dT of the neutral surface at the searched from interface
  real,               intent(in) :: beta_ref    !< dRho/dS of the neutral surface at the searched from interface
  real,               intent(in) :: P_top       !< Pressure at the top interface in the layer to be searched
  real,               intent(in) :: P_bot       !< Pressure at the bottom interface in the layer to be searched
  integer,            intent(in) :: deg         !< Order of the polynomimal used for reconstructions
  real, dimension(:), intent(in) :: ppoly_T     !< Coefficients of the order N polynomial reconstruction of T within
                                                !! the layer to be searched.
  real, dimension(:), intent(in) :: ppoly_S     !< Coefficients of the order N polynomial reconstruction of T within
                                                !! the layer to be searched.
  real,               intent(in) :: x0          !< Nondimensional position within the layer where the neutral
                                                !! surface connects. If interpolate_for_nondim_position was
                                                !! previously called, this would be based on linear profile of dRho
  real,               intent(in) :: drho_top, drho_bot, min_bound
  real,               intent(in) :: ref_pres    !< Optionally use a different reference pressure other than local
  type(EOS_type),     pointer    :: EOS         !< Equation of state structure
  logical, optional,  intent(in) :: force_brent !< Forces the use of Brent's method instead of Newton's method to find
                                                !! position of neutral surface

  ! Local variables
  integer :: form_of_EOS
  integer :: iter
  logical :: do_newton, do_brent

  real :: delta_rho, d_delta_rho_dP ! Terms for the Newton iteration
  real :: P_int, P_min, P_ref ! Interpolated pressure
  real :: delta_rho_init, delta_rho_final
  real :: T, S, alpha, beta, alpha_avg, beta_avg
  ! Newton's Method with variables
  real :: dT_dP, dS_dP, delta_T, delta_S, delta_P
  real :: dbeta_dS, dbeta_dT, dalpha_dT, dalpha_dS, dbeta_dP, dalpha_dP
  real :: a, b, c, b_last
  ! Extra Brent's Method variables
  real :: d, e, f, fa, fb, fc, m, p, q, r, s0, sa, sb, tol, machep

  real :: P_last
  logical :: debug = .false.
  if (ref_pres>=0.) P_ref = ref_pres
  delta_P = P_bot-P_top
  refine_nondim_position = min_bound

  call extract_member_EOS(EOS, form_of_EOS = form_of_EOS)
  do_newton = (form_of_EOS == EOS_LINEAR) .or. (form_of_EOS == EOS_TEOS10) .or. (form_of_EOS == EOS_WRIGHT)
  do_brent = .not. do_newton
  if (present(force_brent)) then
    do_newton = .not. force_brent
    do_brent = force_brent
  endif

  ! Calculate the initial values
  call calc_delta_rho(deg, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S, min_bound, &
                      ref_pres, EOS, delta_rho, P_int, T, S, alpha_avg, beta_avg, delta_T, delta_S)
  delta_rho_init = delta_rho
  if ( ABS(delta_rho_init) < tolerance ) then
    refine_nondim_position = min_bound
    return
  endif

  if ( delta_rho_init > 0.) then
    refine_nondim_position = 1.
    return
  endif

  if (debug) then
    write (*,*) "------"
    write (*,*) "Starting x0, delta_rho: ", min_bound, delta_rho
  endif

  ! For now only linear, Wright, and TEOS-10 equations of state have functions providing second derivatives and
  ! thus can use Newton's method for the equation of state
  if (do_newton) then
    refine_nondim_position = min_bound
    ! Set lower bound of pressure
    P_min = P_top*(1.-min_bound) + P_bot*(min_bound)
    fa = delta_rho_init ; a = min_bound
    fb = delta_rho_init ; b = min_bound
    fc = drho_bot       ; c = 1.
    ! Iterate over Newton's method for the function: x0 = x0 - delta_rho/d_delta_rho_dP
    do iter = 1, max_iter
      P_int = P_top*(1. - b) + P_bot*b
      ! Evaluate total derivative of delta_rho
      if (ref_pres<0.) P_ref = P_int
      call calculate_density_second_derivs( T, S, P_ref, dbeta_dS, dbeta_dT, dalpha_dT, dbeta_dP, dalpha_dP, EOS )
      ! In the case of a constant reference pressure, no dependence on neutral direction with pressure
      if (ref_pres>=0.) then
        dalpha_dP = 0. ; dbeta_dP = 0.
      endif
      dalpha_dS = dbeta_dT ! Cross derivatives are identicial
      ! By chain rule dT_dP= (dT_dz)*(dz/dP) = dT_dz / (Pbot-Ptop)
      dT_dP = first_derivative_polynomial( ppoly_T, deg+1, refine_nondim_position ) / delta_P
      dS_dP = first_derivative_polynomial( ppoly_S, deg+1, refine_nondim_position ) / delta_P
      ! Total derivative of d_delta_rho wrt P
      d_delta_rho_dP = 0.5*( delta_S*(dS_dP*dbeta_dS + dT_dP*dbeta_dT + dbeta_dP) +     &
                             ( delta_T*(dS_dP*dalpha_dS + dT_dP*dalpha_dT + dalpha_dP))) + &
                             dS_dP*beta_avg + dT_dP*alpha_avg
      ! This probably won't happen, but if it does nudge the value a little for the next iteration
      if (d_delta_rho_dP == 0.) then
        b = b + 2*EPSILON(b)*b
      else
        ! Newton step update
        P_int = P_int - (fb / d_delta_rho_dP)
        ! This line is equivalent to the next
        ! refine_nondim_position = (P_top-P_int)/(P_top-P_bot)
        b_last = b
        b = (P_int-P_top)/delta_P
        ! Test to see if it fell out of the bracketing interval. If so, take a bisection step
        if (b < a .or. b > c) then
          b = 0.5*(a + c)
        endif
      endif
      call calc_delta_rho(deg, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S,   &
                          b, ref_pres, EOS, fb, P_int, T, S, alpha_avg, beta_avg, delta_T, delta_S)
      if (debug) print *, "Iteration, b, fb: ", iter, b, fb
      if (ABS(fb) <= tolerance .or. ABS(b-b_last) <= tolerance ) then
        refine_nondim_position = P_int/delta_P
        exit
      endif

      ! Update the bracket
      if (SIGN(1.,fa)*SIGN(1.,fb)<0.) then
        c = b
        fc = delta_rho
      else
        a = b
        fa = delta_rho
      endif
    enddo
    refine_nondim_position = b
    delta_rho = fb
  endif

  ! Do Brent if analytic second derivatives don't exist
  if (do_brent) then
    sa = max(refine_nondim_position,min_bound) ; fa = delta_rho
    sb = 1. ; fb = drho_bot
    c = sa ; fc = fa ; e = sb - sa; d = e


    ! This is from https://people.sc.fsu.edu/~jburkardt/f_src/brent/brent.f90
    do iter = 1,max_iter
      if ( abs ( fc ) < abs ( fb ) ) then
        sa = sb
        sb = c
        c = sa
        fa = fb
        fb = fc
        fc = fa
      end if
      tol = 2. * machep * abs ( sb ) + tolerance
      m = 0.5 * ( c - sb )
      if ( abs ( m ) <= tol .or. fb == 0. ) then
        exit
      endif
      if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then
        e = m
        d = e
      else
        s0 = fb / fa
        if ( sa == c ) then
          p = 2. * m * s0
          q = 1. - s0
        else
          q = fa / fc
          r = fb / fc
          p = s0 * ( 2. * m * q * ( q - r ) - ( sb - sa ) * ( r - 1. ) )
          q = ( q - 1. ) * ( r - 1. ) * ( s0 - 1. )
        end if
        if ( 0. < p ) then
          q = - q
        else
          p = - p
        end if
        s0 = e
        e = d
        if ( 2. * p < 3. * m * q - abs ( tol * q ) .and. &
          p < abs ( 0.5 * s0 * q ) ) then
          d = p / q
        else
          e = m
          d = e
        end if
      end if
      sa = sb
      fa = fb
      if ( tol < abs ( d ) ) then
        sb = sb + d
      else if ( 0. < m ) then
        sb = sb + tol
      else
        sb = sb - tol
      end if
      call calc_delta_rho(deg, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S, &
                        sb, ref_pres, EOS, fb)
      if ( ( 0. < fb .and. 0. < fc ) .or. &
           ( fb <= 0. .and. fc <= 0. ) ) then
        c = sa
        fc = fa
        e = sb - sa
        d = e
      end if
    enddo
    ! Modified from original to ensure that the minimum is found
    fa = ABS(fa) ; fb = ABS(fb) ; fc = ABS(fc)
    delta_rho = MIN(fa, fb, fc)

    if (fb==delta_rho) then
      refine_nondim_position = max(sb,min_bound)
    elseif (fa==delta_rho) then
      refine_nondim_position = max(sa,min_bound)
    elseif (fc==delta_rho) then
      refine_nondim_position = max(c, min_bound)
    endif
  endif

  ! Make sure that the result is bounded between 0 and 1
  if (refine_nondim_position>1.) then
    if (debug) then
      write (*,*) "T, T Poly Coeffs: ", T, ppoly_T
      write (*,*) "S, S Poly Coeffs: ", S, ppoly_S
      write (*,*) "T_ref, alpha_ref: ", T_ref, alpha_ref
      write (*,*) "S_ref, beta_ref : ", S_ref, beta_ref
      write (*,*) "P, dT_dP, dS_dP:", P_int, dT_dP, dS_dP
      write (*,*) "x0: ", x0
      write (*,*) "refine_nondim_position: ", refine_nondim_position
    endif
    call MOM_error(WARNING, "refine_nondim_position>1.")
    refine_nondim_position = MAX(x0,min_bound)
  endif

  if (refine_nondim_position<min_bound) then
    if (debug) then
      write (*,*) "T, T Poly Coeffs: ", T, ppoly_T
      write (*,*) "S, S Poly Coeffs: ", S, ppoly_S
      write (*,*) "T_ref, alpha_ref: ", T_ref, alpha_ref
      write (*,*) "S_ref, beta_ref : ", S_ref, beta_ref
      write (*,*) "dT_dP, dS_dP:", dT_dP, dS_dP
      write (*,*) "x0: ", x0
      write (*,*) "refine_nondim_position: ", refine_nondim_position
    endif
    call MOM_error(WARNING, "refine_nondim_position<min_bound.")
    refine_nondim_position = MAX(x0,min_bound)
  endif

  if (debug) then
    call calc_delta_rho(deg, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S, &
                        refine_nondim_position, ref_pres, EOS, delta_rho)
    write (*,*) "End delta_rho: ", delta_rho
    write (*,*) "x0, delta_x: ", x0, refine_nondim_position-x0
    write (*,*) "refine_nondim_position: ", refine_nondim_position
    write (*,*) "Iterations: ", iter
    write (*,*) "******"
  endif

end function refine_nondim_position

!> Returns .true. if the endpoints of neutral surface do not have the same density (within a specified tolerance)
logical function check_neutral_positions(deg, EOS, x_l, T_poly_l, S_poly_l, P_l, x_r, T_poly_r, S_poly_r, P_r, tolerance, ref_pres)
  integer                 :: deg       !< Degree of polynomial
  type(EOS_type), pointer :: EOS
  real                    :: x_l       !< Nondim position within layer (left)
  real, dimension(deg+1)  :: T_poly_l  !< Coefficients of polynomial reconstructions of T (left)
  real, dimension(deg+1)  :: S_poly_l  !< Coefficients of polynomial reconstructions of S (left)
  real, dimension(2)      :: P_l       !< Pressure at top and bottom of layer (left)
  real                    :: x_r       !< Nondim position within layer (left)
  real, dimension(deg+1)  :: T_poly_r  !< Coefficients of polynomial reconstructions of T (right)
  real, dimension(deg+1)  :: S_poly_r  !< Coefficients of polynomial reconstructions of S (right)
  real, dimension(2)      :: P_r       !< Pressure at top and bottom of layer (right)
  real                    :: tolerance !< How close to the difference in density should be
  real, optional          :: ref_pres  !< reference pressure if not usign local pressure

  real :: delta_rho
  real :: Pl, Tl, Sl, alpha_l, beta_l
  real :: Pr, Tr, Sr, alpha_r, beta_r

  Tl = evaluation_polynomial( T_poly_l, deg+1, x_l )
  Tr = evaluation_polynomial( T_poly_r, deg+1, x_r )
  Sl = evaluation_polynomial( S_poly_l, deg+1, x_l )
  Sr = evaluation_polynomial( S_poly_r, deg+1, x_r )

  if (ref_pres>0.) then
    call calculate_density_derivs( Tl, Sl, ref_pres, alpha_l, beta_l, EOS )
    call calculate_density_derivs( Tr, Sr, ref_pres, alpha_r, beta_r, EOS )
  else
    Pl = (1. - x_l)*P_l(1) + x_l*P_l(2)
    Pr = (1. - x_r)*P_r(1) + x_l*P_r(2)
    call calculate_density_derivs( Tl, Sl, Pl, alpha_l, beta_l, EOS )
    call calculate_density_derivs( Tr, Sr, Pr, alpha_r, beta_r, EOS )
  endif

  delta_rho = 0.5*( (alpha_l+alpha_r)*(Tl-Tr) + (beta_l+beta_r)*(Sl-Sr) )
  check_neutral_positions = ABS(delta_rho)>tolerance

  if (check_neutral_positions) then
    write (*,*) "Density difference of", delta_rho
  endif

end function check_neutral_positions
!> Do a compensated sum to account for roundoff level
subroutine kahan_sum(sum, summand, c)
  real, intent(inout) :: sum      !< Running sum
  real, intent(in   ) :: summand  !< Term to be added
  real ,intent(inout) :: c        !< Keep track of roundoff
  real :: y, t
  y = summand - c
  t = sum + y
  c = (t-sum) - y
  sum = t

end subroutine kahan_sum

end module MOM_neutral_diffusion_aux
