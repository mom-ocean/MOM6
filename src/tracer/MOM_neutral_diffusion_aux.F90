!> A column-wise toolbox for implementing neutral diffusion
module MOM_neutral_diffusion_aux

use MOM_EOS,                   only : EOS_type, extract_member_EOS, EOS_LINEAR, EOS_TEOS10, EOS_WRIGHT
use MOM_EOS,                   only : calculate_density_derivs, calculate_density_second_derivs
use MOM_error_handler,         only : MOM_error, FATAL, WARNING
use polynomial_functions,      only : evaluation_polynomial, first_derivative_polynomial

! This file is part of MOM6. See LICENSE.md for the license.
implicit none ; private

public set_ndiff_aux_params
public mark_unstable_cells
public increment_interface
public calc_drho
public drho_at_pos
public search_other_column
public interpolate_for_nondim_position
public refine_nondim_position
public check_neutral_positions
public kahan_sum

!> The control structure for this module
type, public :: ndiff_aux_CS_type ; private
  integer :: nterm        !< Number of terms in polynomial (deg+1)
  integer :: max_iter     !< Maximum number of iterations
  real :: drho_tol        !< Tolerance criterion for difference in density (kg/m3)
  real :: xtol            !< Criterion for how much position changes (nondim)
  real :: ref_pres        !< Determines whether a constant reference pressure is used everywhere or locally referenced
                          !< density is done. ref_pres <-1 is the latter, ref_pres >= 0. otherwise
  logical :: force_brent = .false. !< Use Brent's method instead of Newton even when second derivatives are available
  logical :: debug        !< If true, write verbose debugging messages and checksusm
  type(EOS_type), pointer :: EOS !< Pointer to equation of state used in the model
end type ndiff_aux_CS_type

contains

!> Initialize the parameters used to iteratively find the neutral direction
subroutine set_ndiff_aux_params( CS, deg, max_iter, drho_tol, xtol, ref_pres, force_brent, EOS, debug)
  type(ndiff_aux_CS_type),  intent(inout) :: CS          !< Control structure for refine_pos
  integer,        optional, intent(in   ) :: deg         !< Degree of polynommial used in reconstruction
  integer,        optional, intent(in   ) :: max_iter    !< Maximum number of iterations
  real,           optional, intent(in   ) :: drho_tol    !< Tolerance for function convergence
  real,           optional, intent(in   ) :: xtol        !< Tolerance for change in position
  real,           optional, intent(in   ) :: ref_pres    !< Reference pressure to use
  logical,        optional, intent(in   ) :: force_brent !< Force Brent method for linear, TEOS-10, and WRIGHT
  logical,        optional, intent(in   ) :: debug       !< If true, print output use to help debug neutral diffusion
  type(EOS_type), target, optional, intent(in   ) :: EOS !< Equation of state

  if (present( deg         )) CS%nterm       =  deg + 1
  if (present( max_iter    )) CS%max_iter    =  max_iter
  if (present( drho_tol    )) CS%drho_tol    =  drho_tol
  if (present( xtol        )) CS%xtol        =  xtol
  if (present( ref_pres    )) CS%ref_pres    =  ref_pres
  if (present( force_brent )) CS%force_brent =  force_brent
  if (present( EOS         )) CS%EOS         => EOS
  if (present( debug       )) CS%debug       =  debug

end subroutine set_ndiff_aux_params

!> Given the reconsturcitons of dRdT, dRdS, T, S mark the cells which are stably stratified parts of the water column
!! For an layer to be unstable the top interface must be denser than the bottom or the bottom interface of the layer
subroutine mark_unstable_cells(nk, dRdT, dRdS,T, S, stable_cell, ns)
  integer,                intent(in)    :: nk          !< Number of levels in a column
  real, dimension(nk,2),  intent(in)    :: dRdT        !< drho/dT (kg/m3/degC) at interfaces
  real, dimension(nk,2),  intent(in)    :: dRdS        !< drho/dS (kg/m3/ppt) at interfaces
  real, dimension(nk,2),  intent(in)    :: T           !< drho/dS (kg/m3/ppt) at interfaces
  real, dimension(nk,2),  intent(in)    :: S           !< drho/dS (kg/m3/ppt) at interfaces
  logical, dimension(nk), intent(  out) :: stable_cell !< True if this cell is unstably stratified
  integer,                intent(  out) :: ns          !< Number of neutral surfaces in unmasked part of the column

  integer :: k, first_stable, prev_stable
  real :: delta_rho

  ! First check to make sure that density profile between the two interfaces of the cell are stable
  ! Note that we neglect a factor of 0.5 because we only care about the sign of delta_rho not magnitude
  do k = 1,nk
    ! Compare density of bottom interface to top interface, should be positive (or zero) if stable
    delta_rho = (dRdT(k,2) + dRdT(k,1))*(T(k,2) - T(k,1)) + (dRdS(k,2) + dRdS(k,1))*(S(k,2) - S(k,1))
    stable_cell(k) = delta_rho >= 0.
  enddo

  first_stable = 1
  ! Check to see that bottom interface of upper cell is lighter than the upper interface of the lower cell
  do k=1,nk
    if (stable_cell(k)) then
      first_stable = k
      exit
    endif
  enddo
  prev_stable = first_stable

  ! Start either with the first stable cell or the layer just below the surface
  do k = prev_stable+1, nk
    ! Don't do anything if the cell has already been marked as unstable
    if (.not. stable_cell(k)) cycle
    ! Otherwise, we need to check to see if this cell's upper interface is denser than the previous stable_cell
    ! Compare top interface of lower cell to bottom interface of upper cell, positive or zero if bottom cell is stable
    delta_rho = (dRdT(k,1) + dRdT(prev_stable,2))*(T(k,1) - T(prev_stable,2)) + &
                (dRdS(k,1) + dRdS(prev_stable,2))*(S(k,1) - S(prev_stable,2))
    stable_cell(k) = delta_rho >= 0.
    ! If the lower cell is marked as stable, then it should be the next reference cell
    if (stable_cell(k)) prev_stable = k
  enddo

  ! Number of interfaces is the 2 times number of stable cells in the water column
  ns = 0
  do k = 1,nk
    if (stable_cell(k)) ns = ns + 2
  enddo

end subroutine mark_unstable_cells

!> Increments the interface which was just connected and also set flags if the bottom is reached
subroutine increment_interface(nk, kl, ki, stable, reached_bottom, searching_this_column, searching_other_column)
  integer, intent(in   )                :: nk                     !< Number of vertical levels
  integer, intent(inout)                :: kl                     !< Current layer (potentially updated)
  integer, intent(inout)                :: ki                     !< Current interface
  logical, dimension(nk), intent(in   ) :: stable                 !< True if the cell is stably stratified
  logical, intent(inout)                :: reached_bottom         !< Updated when kl == nk and ki == 2
  logical, intent(inout)                :: searching_this_column  !< Updated when kl == nk and ki == 2
  logical, intent(inout)                :: searching_other_column !< Updated when kl == nk and ki == 2
  integer :: k

  if (ki == 1) then
    ki = 2
  elseif ((ki == 2) .and. (kl < nk) ) then
    do k = kl+1,nk
      if (stable(kl)) then
        kl = k
        ki = 1
        exit
      endif
      ! If we did not find another stable cell, then the current cell is essentially the bottom
      ki = 2
      reached_bottom = .true.
      searching_this_column = .true.
      searching_other_column = .false.
    enddo
  elseif ((kl == nk) .and. (ki==2)) then
    reached_bottom = .true.
    searching_this_column = .true.
    searching_other_column = .false.
  else
    call MOM_error(FATAL,"Unanticipated eventuality in increment_interface")
  endif
end subroutine increment_interface

!> Calculates difference in density at two points (rho1-rho2) with known density derivatives, T, and S
real function calc_drho(T1, S1, dRdT1, dRdS1, T2, S2, dRdT2, dRdS2)
  real, intent(in   ) :: T1     !< Temperature at point 1
  real, intent(in   ) :: S1     !< Salinity at point 1
  real, intent(in   ) :: dRdT1  !< dRhodT at point 1
  real, intent(in   ) :: dRdS1  !< dRhodS at point 1
  real, intent(in   ) :: T2     !< Temperature at point 2
  real, intent(in   ) :: S2     !< Salinity at point 2
  real, intent(in   ) :: dRdT2  !< dRhodT at point 2
  real, intent(in   ) :: dRdS2  !< dRhodS at point

  calc_drho = 0.5*( (dRdT1+dRdT2)*(T1-T2) + (dRdS1+dRdS2)*(S1-S2) )
end function calc_drho

!> Calculate the difference in neutral density between a reference T, S, alpha, and beta
!! at a point on the polynomial reconstructions of T, S
subroutine drho_at_pos(CS, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S, x0, &
                       delta_rho, P_out, T_out, S_out, alpha_avg_out, beta_avg_out, delta_T_out, delta_S_out)
  type(ndiff_aux_CS_type), intent(in) :: CS           !< Control structure with parameters for this module
  real,                   intent(in)  :: T_ref     !< Temperature at reference surface
  real,                   intent(in)  :: S_ref     !< Salinity at reference surface
  real,                   intent(in)  :: alpha_ref !< dRho/dT at reference surface
  real,                   intent(in)  :: beta_ref  !< dRho/dS at reference surface
  real,                   intent(in)  :: P_top     !< Pressure (Pa) at top interface of layer to be searched
  real,                   intent(in)  :: P_bot     !< Pressure (Pa) at bottom interface
  real, dimension(CS%nterm), intent(in)  :: ppoly_T   !< Coefficients of T reconstruction
  real, dimension(CS%nterm), intent(in)  :: ppoly_S   !< Coefficients of S reconstruciton
  real,                   intent(in)  :: x0        !< Nondimensional position to evaluate
  real,                   intent(out) :: delta_rho !< The density difference from a reference value
  real,         optional, intent(out) :: P_out         !< Pressure at point x0
  real,         optional, intent(out) :: T_out         !< Temperature at point x0
  real,         optional, intent(out) :: S_out         !< Salinity at point x0
  real,         optional, intent(out) :: alpha_avg_out !< Average of alpha between reference and x0
  real,         optional, intent(out) :: beta_avg_out  !< Average of beta between reference and x0
  real,         optional, intent(out) :: delta_T_out   !< Difference in temperature between reference and x0
  real,         optional, intent(out) :: delta_S_out   !< Difference in salinity between reference and x0

  real :: alpha, beta, alpha_avg, beta_avg, P_int, T, S, delta_T, delta_S

  P_int = (1. - x0)*P_top + x0*P_bot
  T = evaluation_polynomial( ppoly_T, CS%nterm, x0 )
  S = evaluation_polynomial( ppoly_S, CS%nterm, x0 )
  ! Interpolated pressure if using locally referenced neutral density
  if (CS%ref_pres<0.) then
    call calculate_density_derivs( T, S, P_int, alpha, beta, CS%EOS )
  else
  ! Constant reference pressure (isopycnal)
    call calculate_density_derivs( T, S, CS%ref_pres, alpha, beta, CS%EOS )
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

end subroutine drho_at_pos

!> Searches the "other" (searched) column for the position of the neutral surface
subroutine search_other_column(dRhoTop, dRhoBot, Ptop, Pbot, lastP, lastK, kl, kl_0, ki, &
                               top_connected, bot_connected, out_P, out_K, search_layer)
  real,                  intent(in   ) :: dRhoTop        !< Density difference across top interface
  real,                  intent(in   ) :: dRhoBot        !< Density difference across top interface
  real,                  intent(in   ) :: Ptop           !< Pressure at top interface
  real,                  intent(in   ) :: Pbot           !< Pressure at bottom interface
  real,                  intent(in   ) :: lastP          !< Last position connected in the searched column
  integer,               intent(in   ) :: lastK          !< Last layer connected in the searched column
  integer,               intent(in   ) :: kl             !< Layer in the searched column
  integer,               intent(in   ) :: kl_0           !< Layer in the searched column
  integer,               intent(in   ) :: ki             !< Interface of the searched column
  logical, dimension(:), intent(inout) :: top_connected  !< True if the top interface was pointed to
  logical, dimension(:), intent(inout) :: bot_connected  !< True if the top interface was pointed to
  real,                  intent(  out) :: out_P          !< Position within searched column
  integer,               intent(  out) :: out_K          !< Layer within searched column
  logical,               intent(  out) :: search_layer   !< Neutral surface within cell

  search_layer = .false.
  if (kl > kl_0) then ! Away from top cell
    if (kl == lastK) then ! Searching in the same layer
      if (dRhoTop > 0.) then
        if (lastK == kl) then
          out_P = lastP
        else
          out_P = 0.
        endif
        out_K = kl
!        out_P = max(0.,lastP) ; out_K = kl
      elseif ( dRhoTop == dRhoBot ) then
        if (top_connected(kl)) then
          out_P = 1. ; out_K = kl
        else
          out_P = max(0.,lastP) ; out_K = kl
        endif
      elseif (dRhoTop >= dRhoBot) then
        out_P = 1. ; out_K = kl
      elseif (dRhoTop < 0. .and. dRhoBot < 0.)then
        out_P = 1. ; out_K = kl
      else
        out_K = kl
        out_P = max(interpolate_for_nondim_position( dRhoTop, Ptop, dRhoBot, Pbot ),lastP)
        search_layer = .true.
      endif
    else ! Searching across the interface
      if (.not. bot_connected(kl-1) ) then
        out_K = kl-1
        out_P = 1.
      else
        out_K = kl
        out_P = 0.
      endif
    endif
  else ! At the top cell
    if (ki == 1) then
      out_P = 0. ; out_K = kl
    elseif (dRhoTop > 0.) then
      if (lastK == kl) then
        out_P = lastP
      else
        out_P = 0.
      endif
      out_K = kl
!      out_P = max(0.,lastP) ; out_K = kl
    elseif ( dRhoTop == dRhoBot ) then
      if (top_connected(kl)) then
        out_P = 1. ; out_K = kl
      else
        out_P = max(0.,lastP) ; out_K = kl
      endif
    elseif (dRhoTop >= dRhoBot) then
      out_P = 1. ; out_K = kl
    elseif (dRhoTop < 0. .and. dRhoBot < 0.)then
      out_P = 1. ; out_K = kl
    else
      out_K = kl
      out_P = max(interpolate_for_nondim_position( dRhoTop, Ptop, dRhoBot, Pbot ),lastP)
      search_layer = .true.
    endif
  endif

end subroutine search_other_column

!> Returns the non-dimensional position between Pneg and Ppos where the
!! interpolated density difference equals zero.
!! The result is always bounded to be between 0 and 1.
real function interpolate_for_nondim_position(dRhoNeg, Pneg, dRhoPos, Ppos)
  real, intent(in) :: dRhoNeg !< Negative density difference
  real, intent(in) :: Pneg    !< Position of negative density difference
  real, intent(in) :: dRhoPos !< Positive density difference
  real, intent(in) :: Ppos    !< Position of positive density difference

  if (Ppos<Pneg) then
    stop 'interpolate_for_nondim_position: Houston, we have a problem! Ppos<Pneg'
  elseif (dRhoNeg>dRhoPos) then
    write(0,*) 'dRhoNeg, Pneg, dRhoPos, Ppos=',dRhoNeg, Pneg, dRhoPos, Ppos
  elseif (dRhoNeg>dRhoPos) then
    stop 'interpolate_for_nondim_position: Houston, we have a problem! dRhoNeg>dRhoPos'
  endif
  if (Ppos<=Pneg) then ! Handle vanished or inverted layers
    interpolate_for_nondim_position = 0.5
  elseif ( dRhoPos - dRhoNeg > 0. ) then
    interpolate_for_nondim_position = min( 1., max( 0., -dRhoNeg / ( dRhoPos - dRhoNeg ) ) )
  elseif ( dRhoPos - dRhoNeg == 0) then
    if (dRhoNeg>0.) then
      interpolate_for_nondim_position = 0.
    elseif (dRhoNeg<0.) then
      interpolate_for_nondim_position = 1.
    else ! dRhoPos = dRhoNeg = 0
      interpolate_for_nondim_position = 0.5
    endif
  else ! dRhoPos - dRhoNeg < 0
    interpolate_for_nondim_position = 0.5
  endif
  if ( interpolate_for_nondim_position < 0. ) &
    stop 'interpolate_for_nondim_position: Houston, we have a problem! Pint < Pneg'
  if ( interpolate_for_nondim_position > 1. ) &
    stop 'interpolate_for_nondim_position: Houston, we have a problem! Pint > Ppos'
end function interpolate_for_nondim_position

!> Use root-finding methods to find where dRho = 0, based on the equation of state and the polynomial
!! reconstructions of temperature, salinity. Initial guess is based on the zero crossing of based on linear
!! profiles of dRho, T, and S, between the top and bottom interface. If second derivatives of the EOS are available,
!! it starts with a Newton's method. However, Newton's method is not guaranteed to be bracketed, a check is performed
!! to see if it it diverges outside the interval. In that case (or in the case that second derivatives are not
!! available), Brent's method is used following the implementation found at
!! https://people.sc.fsu.edu/~jburkardt/f_src/brent/brent.f90
real function refine_nondim_position(CS, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, &
                                     ppoly_T, ppoly_S, drho_top, drho_bot, min_bound)
  type(ndiff_aux_CS_type),  intent(in) :: CS        !< Control structure with parameters for this module
  real,                     intent(in) :: T_ref     !< Temperature of the neutral surface at the searched from interface
  real,                     intent(in) :: S_ref     !< Salinity of the neutral surface at the searched from interface
  real,                     intent(in) :: alpha_ref !< dRho/dT of the neutral surface at the searched from interface
  real,                     intent(in) :: beta_ref  !< dRho/dS of the neutral surface at the searched from interface
  real,                     intent(in) :: P_top     !< Pressure at the top interface in the layer to be searched
  real,                     intent(in) :: P_bot     !< Pressure at the bottom interface in the layer to be searched
  real, dimension(:),       intent(in) :: ppoly_T   !< Coefficients of the order N polynomial reconstruction of T within
                                                    !! the layer to be searched.
  real, dimension(:),       intent(in) :: ppoly_S   !< Coefficients of the order N polynomial reconstruction of T within
                                                    !! the layer to be searched.
  real,                     intent(in) :: drho_top  !< Delta rho at top interface (or previous position in cell
  real,                     intent(in) :: drho_bot  !< Delta rho at bottom interface
  real,                     intent(in) :: min_bound !< Lower bound of position, also serves as the initial guess

  ! Local variables
  integer :: form_of_EOS
  integer :: iter
  logical :: do_newton, do_brent

  real :: delta_rho, d_delta_rho_dP ! Terms for the Newton iteration
  real :: P_int, P_min, P_ref ! Interpolated pressure
  real :: delta_rho_init, delta_rho_final
  real :: neg_x, neg_fun
  real :: T, S, alpha, beta, alpha_avg, beta_avg
  ! Newton's Method with variables
  real :: dT_dP, dS_dP, delta_T, delta_S, delta_P
  real :: dbeta_dS, dbeta_dT, dalpha_dT, dalpha_dS, dbeta_dP, dalpha_dP
  real :: a, b, c, b_last
  ! Extra Brent's Method variables
  real :: d, e, f, fa, fb, fc, m, p, q, r, s0, sa, sb, tol, machep

  real :: P_last

  machep = EPSILON(machep)
  if (CS%ref_pres>=0.) P_ref = CS%ref_pres
  delta_P = P_bot-P_top
  refine_nondim_position = min_bound

  call extract_member_EOS(CS%EOS, form_of_EOS = form_of_EOS)
  do_newton = (form_of_EOS == EOS_LINEAR) .or. (form_of_EOS == EOS_TEOS10) .or. (form_of_EOS == EOS_WRIGHT)
  do_brent = .not. do_newton
  if (CS%force_brent) then
    do_newton = .not. CS%force_brent
    do_brent = CS%force_brent
  endif

  ! Calculate the initial values
  call drho_at_pos(CS, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S, min_bound, &
                   delta_rho, P_int, T, S, alpha_avg, beta_avg, delta_T, delta_S)
  delta_rho_init = delta_rho
  if ( ABS(delta_rho_init) <= CS%drho_tol ) then
    refine_nondim_position = min_bound
    return
  endif
  if (ABS(drho_bot) <= CS%drho_tol) then
    refine_nondim_position = 1.
    return
  endif

  ! Set the initial values to ensure that the algorithm returns a 'negative' value
  neg_fun = delta_rho
  neg_x = min_bound

  if (CS%debug) then
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
    do iter = 1, CS%max_iter
      P_int = P_top*(1. - b) + P_bot*b
      ! Evaluate total derivative of delta_rho
      if (CS%ref_pres<0.) P_ref = P_int
      call calculate_density_second_derivs( T, S, P_ref, dbeta_dS, dbeta_dT, dalpha_dT, dbeta_dP, dalpha_dP, CS%EOS )
      ! In the case of a constant reference pressure, no dependence on neutral direction with pressure
      if (CS%ref_pres>=0.) then
        dalpha_dP = 0. ; dbeta_dP = 0.
      endif
      dalpha_dS = dbeta_dT ! Cross derivatives are identicial
      ! By chain rule dT_dP= (dT_dz)*(dz/dP) = dT_dz / (Pbot-Ptop)
      dT_dP = first_derivative_polynomial( ppoly_T, CS%nterm, b ) / delta_P
      dS_dP = first_derivative_polynomial( ppoly_S, CS%nterm, b ) / delta_P
      ! Total derivative of d_delta_rho wrt P
      d_delta_rho_dP = 0.5*( delta_S*(dS_dP*dbeta_dS + dT_dP*dbeta_dT + dbeta_dP) +     &
                             ( delta_T*(dS_dP*dalpha_dS + dT_dP*dalpha_dT + dalpha_dP))) + &
                             dS_dP*beta_avg + dT_dP*alpha_avg
      ! This probably won't happen, but if it does take a bisection
      if (d_delta_rho_dP == 0.) then
        b = 0.5*(a+c)
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
      call drho_at_pos(CS, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S,   &
                       b, fb, P_int, T, S, alpha_avg, beta_avg, delta_T, delta_S)
      if (CS%debug) write(*,'(A,I3.3,X,ES23.15,X,ES23.15)') "Iteration, b, fb: ", iter, b, fb

      if (fb < 0. .and. fb > neg_fun) then
        neg_fun = fb
        neg_x = b
      endif

      ! For the logic to find neutral surfaces to work properly, the function needs to converge to zero
      !  or a small negative value
      if ((fb <= 0.) .and. (fb >= -CS%drho_tol)) then
        refine_nondim_position = b
        exit
      endif
      ! Exit if method has stalled out
      if ( ABS(b-b_last)<=CS%xtol ) then
        refine_nondim_position = b
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
  if (delta_rho > 0.) then
    refine_nondim_position = neg_x
    delta_rho = neg_fun
  endif
  ! Do Brent if analytic second derivatives don't exist
  if (do_brent) then
    sa = max(refine_nondim_position,min_bound) ; fa = delta_rho
    sb = 1. ; fb = drho_bot
    c = sa ; fc = fa ; e = sb - sa; d = e


    ! This is from https://people.sc.fsu.edu/~jburkardt/f_src/brent/brent.f90
    do iter = 1,CS%max_iter
      if ( abs ( fc ) < abs ( fb ) ) then
        sa = sb
        sb = c
        c = sa
        fa = fb
        fb = fc
        fc = fa
      endif
      tol = 2. * machep * abs ( sb ) + CS%xtol
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
        endif
        if ( 0. < p ) then
          q = - q
        else
          p = - p
        endif
        s0 = e
        e = d
        if ( 2. * p < 3. * m * q - abs ( tol * q ) .and. &
          p < abs ( 0.5 * s0 * q ) ) then
          d = p / q
        else
          e = m
          d = e
        endif
      endif
      sa = sb
      fa = fb
      if ( tol < abs ( d ) ) then
        sb = sb + d
      elseif ( 0. < m ) then
        sb = sb + tol
      else
        sb = sb - tol
      endif
      call drho_at_pos(CS, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S, &
                       sb, fb)
      if ( ( 0. < fb .and. 0. < fc ) .or. &
           ( fb <= 0. .and. fc <= 0. ) ) then
        c = sa
        fc = fa
        e = sb - sa
        d = e
      endif
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
    if (CS%debug) then
      write (*,*) "T, T Poly Coeffs: ", T, ppoly_T
      write (*,*) "S, S Poly Coeffs: ", S, ppoly_S
      write (*,*) "T_ref, alpha_ref: ", T_ref, alpha_ref
      write (*,*) "S_ref, beta_ref : ", S_ref, beta_ref
      write (*,*) "x0: ", min_bound
      write (*,*) "refine_nondim_position: ", refine_nondim_position
    endif
    call MOM_error(WARNING, "refine_nondim_position>1.")
    refine_nondim_position = 1.
  endif

  if (refine_nondim_position<min_bound) then
    if (CS%debug) then
      write (*,*) "T, T Poly Coeffs: ", T, ppoly_T
      write (*,*) "S, S Poly Coeffs: ", S, ppoly_S
      write (*,*) "T_ref, alpha_ref: ", T_ref, alpha_ref
      write (*,*) "S_ref, beta_ref : ", S_ref, beta_ref
      write (*,*) "x0: ", min_bound
      write (*,*) "refine_nondim_position: ", refine_nondim_position
    endif
    call MOM_error(WARNING, "refine_nondim_position<min_bound.")
    refine_nondim_position = min_bound
  endif

  if (CS%debug) then
    call drho_at_pos(CS, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S, &
                     refine_nondim_position, delta_rho)
    write (*,*) "End delta_rho: ", delta_rho
    write (*,*) "x0, delta_x: ", min_bound, refine_nondim_position-min_bound
    write (*,*) "refine_nondim_position: ", refine_nondim_position
    write (*,*) "Iterations: ", iter
    write (*,*) "T_poly_coeffs: ", ppoly_T
    write (*,*) "S_poly_coeffs: ", ppoly_S
    write (*,*) "******"
  endif

end function refine_nondim_position

subroutine check_neutral_positions(CS, Ptop_l, Pbot_l, Ptop_r, Pbot_r, PoL, PoR, Tl_coeffs, Tr_coeffs, &
                                 Sl_coeffs, Sr_coeffs)
  type(ndiff_aux_CS_type),   intent(in) :: CS        !< Control structure with parameters for this module
  real,                      intent(in) :: Ptop_l    !< Pressure at top interface of left cell
  real,                      intent(in) :: Pbot_l    !< Pressure at bottom interface of left cell
  real,                      intent(in) :: Ptop_r    !< Pressure at top interface of right cell
  real,                      intent(in) :: Pbot_r    !< Pressure at bottom interface of right cell
  real,                      intent(in) :: PoL       !< Nondim position in left cell
  real,                      intent(in) :: PoR       !< Nondim position in right cell
  real, dimension(CS%nterm), intent(in) :: Tl_coeffs !< T polynomial coefficients of left cell
  real, dimension(CS%nterm), intent(in) :: Tr_coeffs !< T polynomial coefficients of right cell
  real, dimension(CS%nterm), intent(in) :: Sl_coeffs !< S polynomial coefficients of left cell
  real, dimension(CS%nterm), intent(in) :: Sr_coeffs !< S polynomial coefficients of right cell
  ! Local variables
  real :: Pl, Pr, Tl, Tr, Sl, Sr
  real :: alpha_l, alpha_r, beta_l, beta_r
  real :: delta_rho

  Pl = (1. - PoL)*Ptop_l + PoL*Pbot_l
  Pr = (1. - PoR)*Ptop_r + PoR*Pbot_r
  Tl = evaluation_polynomial( Tl_coeffs, CS%nterm, PoL )
  Sl = evaluation_polynomial( Sl_coeffs, CS%nterm, PoL )
  Tr = evaluation_polynomial( Tr_coeffs, CS%nterm, PoR )
  Sr = evaluation_polynomial( Sr_coeffs, CS%nterm, PoR )

  if (CS%ref_pres<0.) then
    call calculate_density_derivs( Tl, Sl, Pl, alpha_l, beta_l, CS%EOS )
    call calculate_density_derivs( Tr, Sr, Pr, alpha_r, beta_r, CS%EOS )
  else
    call calculate_density_derivs( Tl, Sl, CS%ref_pres, alpha_l, beta_l, CS%EOS )
    call calculate_density_derivs( Tr, Sr, CS%ref_pres, alpha_r, beta_r, CS%EOS )
  endif

  delta_rho = 0.5*( (alpha_l + alpha_r)*(Tl - Tr) + (beta_l + beta_r)*(Sl - Sr) )
  write(*,'(A,ES23.15)') "Delta-rho: ", delta_rho

end subroutine check_neutral_positions

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
