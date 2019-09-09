!> Calculate and apply diffusive fluxes as a parameterization of lateral mixing (non-neutral) by
!! mesoscale eddies near the top and bottom boundary layers of the ocean.
module MOM_boundary_lateral_mixing

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,             only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_diag_mediator,         only : post_data, register_diag_field
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_file_parser,           only : openParameterBlock, closeParameterBlock
use MOM_grid,                  only : ocean_grid_type
use MOM_remapping,             only : remapping_CS, initialize_remapping
use MOM_remapping,             only : extract_member_remapping_CS, build_reconstructions_1d
use MOM_remapping,             only : average_value_ppoly, remappingSchemesDoc, remappingDefaultScheme
use MOM_tracer_registry,       only : tracer_registry_type, tracer_type
use MOM_verticalGrid,          only : verticalGrid_type
use polynomial_functions,      only : evaluation_polynomial, first_derivative_polynomial

implicit none ; private

public near_boundary_unit_tests

! Private parameters to avoid doing string comparisons for bottom or top boundary layer
integer, parameter :: SURFACE = -1 !< Set a value that corresponds to the surface bopundary
integer, parameter :: BOTTOM  = 1  !< Set a value that corresponds to the bottom boundary

#include <MOM_memory.h>
contains

!> Driver routine for calculating lateral diffusive fluxes near the top and bottom boundaries. Two different methods
!! Method 1: Calculate fluxes from bulk layer integrated quantities
subroutine boundary_lateral_mixing()


end subroutine

!< Calculate bulk layer value of a scalar quantity as the thickness weighted average
real function bulk_average(h, hBLT, phi)
  real, dimension(:), intent(in) :: h    !< Layer thicknesses [m]
  real              , intent(in) :: hBLT !< Depth of the mixing layer [m]
  real, dimension(:), intent(in) :: phi  !< Scalar quantity
  ! Local variables
  integer :: nk ! Number of layers
  real    :: htot ! Running sum of the thicknesses (top to bottom)
  integer :: k

  ! if ( len(h) .ne. len(phi ) call MOM_error(FATAL,"boundary_mixing: tracer and thicknesses of different size")
  nk = SIZE(h)

  htot = 0.
  bulk_average = 0.
  do k = 1,nk
    bulk_average = bulk_average + phi(k)*h(k)
    htot = htot + h(k)
  enddo

  if (htot > 0.) then
      bulk_average = bulk_average / hBLT
  else
      call MOM_error(FATAL, "Column thickness is 0.")
  endif

end function bulk_average

!> Calculate the harmonic mean of two quantities
real function harmonic_mean(h1,h2)
  real :: h1 !< Scalar quantity
  real :: h2 !< Scalar quantity

  harmonic_mean = 2.*(h1*h2)/(h1+h2)
end function harmonic_mean

!> Find the k-index range corresponding to the layers that are within the boundary-layer region
subroutine boundary_k_range(boundary, nk, h, hbl, k_top, zeta_top, k_bot, zeta_bot)
  integer,             intent(in   ) :: boundary !< SURFACE or BOTTOM                       [nondim]
  integer,             intent(in   ) :: nk       !< Number of layers                        [nondim]
  real, dimension(nk), intent(in   ) :: h        !< Layer thicknesses of the coluymn        [m]
  real,                intent(in   ) :: hbl      !< Thickness of the boundary layer         [m]
                                                 !! If surface, with respect to zbl_ref = 0.
                                                 !! If bottom, with respect to zbl_ref = SUM(h)
  integer,             intent(  out) :: k_top    !< Index of the first layer within the boundary
  real,                intent(  out) :: zeta_top !< Distance from the top of a layer to the intersection of the
                                                 !! top extent of the boundary layer (0 at top, 1 at bottom)  [nondim]
  integer,             intent(  out) :: k_bot    !< Index of the last layer within the boundary
  real,                intent(  out) :: zeta_bot !< Distance of the lower layer to the boundary layer depth
                                                 !! (0 at top, 1 at bottom)  [nondim]
  ! Local variables
  real :: htot
  integer :: k
  ! Surface boundary layer
  if ( boundary == SURFACE ) then
    k_top = 1
    zeta_top = 0.
    htot = 0.
    do k=1,nk
      htot = htot + h(k)
      if ( htot >= hbl) then
        k_bot = k
        zeta_bot = 1 - (htot - hbl)/h(k)
        return
      endif
    enddo
  ! Bottom boundary layer
  elseif ( boundary == BOTTOM ) then
    k_bot = nk
    zeta_bot = 1.
    htot = 0.
    do k=nk,1,-1
      htot = htot + h(k)
      if (htot >= hbl) then
        k_top = k
        zeta_top = 1 - (htot - hbl)/h(k)
        return
      endif
    enddo
  else
    call MOM_error(FATAL,"Houston, we've had a problem in boundary_k_range")
  endif

end subroutine boundary_k_range

!> Calculate the near-boundary diffusive fluxes calculated from a 'bulk model'
subroutine layer_fluxes_bulk_method(boundary, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R, phi_pp_L, phi_pp_R, &
                                    khtr_u, F_layer)
  integer,                   intent(in   )       :: boundary !< Which boundary layer SURFACE or BOTTOM  [nondim]
  integer,                   intent(in   )       :: nk       !< Number of layers                        [nondim]
  integer,                   intent(in   )       :: deg      !< order of the polynomial reconstruction  [nondim]
  real, dimension(nk),       intent(in   )       :: h_L      !< Layer thickness (left)                  [m]
  real, dimension(nk),       intent(in   )       :: h_R      !< Layer thickness (right)                 [m]
  real,                      intent(in   )       :: hbl_L    !< Thickness of the boundary boundary
                                                                       !! layer (left)                            [m]
  real,                      intent(in   )       :: hbl_R    !< Thickness of the boundary boundary
                                                             !! layer (left)                            [m]
  real, dimension(nk),       intent(in   )       :: phi_L    !< Tracer values (left)                    [ nondim m^-3 ]
  real, dimension(nk),       intent(in   )       :: phi_R    !< Tracer values (right)                   [ nondim m^-3 ]
  real, dimension(nk,deg+1), intent(in   )       :: phi_pp_L !< Tracer reconstruction (left)            [ nondim m^-3 ]
  real, dimension(nk,deg+1), intent(in   )       :: phi_pp_R !< Tracer reconstruction (right)           [ nondim m^-3 ]
  real, dimension(nk),       intent(in   )       :: khtr_u   !< Horizontal diffusivities at U-point     [m^2 s^-1]
  real, dimension(nk),       intent(  out)       :: F_layer  !< Layerwise diffusive flux at U-point     [trunit s^-1]
  ! Local variables
  real                :: F_bulk               ! Total diffusive flux across the U point           [trunit s^-1]
  real, dimension(nk) :: h_means              ! Calculate the layer-wise harmonic means           [m]
  real, dimension(nk) :: h_u                  ! Thickness at the u-point                          [m]
  real                :: hbl_u                ! Boundary layer Thickness at the u-point           [m]
  real                :: khtr_avg             ! Thickness-weighted diffusivity at the u-point     [m^2 s^-1]
  real                :: heff                 ! Harmonic mean of layer thicknesses                [m]
  real                :: inv_heff             ! Inverse of the harmonic mean of layer thicknesses [m^[-1]
  real                :: phi_L_avg, phi_R_avg ! Bulk, thickness-weighted tracer averages (left and right column)
                                              !                                                   [trunit m^-3 ]
  real    :: htot                 ! Total column thickness [m]
  integer :: k
  integer :: k_top_L, k_bot_L
  integer :: k_top_R, k_bot_R

  ! Calculate bulk averages of various quantities
  phi_L_avg = bulk_average(h_L, hbl_L, phi_L)
  phi_R_avg = bulk_average(h_R, hbl_R, phi_R)
  do k=1,nk
    h_u(k) = 0.5 * (h_L(k) + h_R(k))
  enddo
  hbl_u = 0.5*(hbl_L + hbl_R)
  khtr_avg  = bulk_average(h_u, hbl_u, khtr_u)

  ! Calculate the 'bulk' diffusive flux from the bulk averaged quantities
  heff = harmonic_mean(hbl_L, hbl_R)
  F_bulk = -(khtr_avg * heff) * (phi_R_avg - phi_L_avg)

  ! Calculate the layerwise sum of the vertical effective thickness. This is different than the heff calculated
  ! above, but is used as a way to decompose decompose the fluxes onto the individual layers
  do k=1,nk
    h_means(k) = harmonic_mean(h_L(k),h_R(k))
  enddo
  inv_heff = 1./SUM(h_means)
  do k=1,nk
    F_layer(k) = F_bulk * (h_means(k)*inv_heff)
  enddo

end subroutine layer_fluxes_bulk_method

!> Unit tests for near-boundary horizontal mixing
logical function near_boundary_unit_tests( verbose )
  logical,               intent(in) :: verbose !< If true, output additional information for debugging unit tests

  ! Local variables
  integer, parameter    :: nk = 2               ! Number of layers
  integer, parameter    :: deg = 1              ! Degree of reconstruction (linear here)
  real, dimension(nk)   :: phi_L, phi_R         ! Tracer values (left and right column)             [ nondim m^-3 ]
  real, dimension(nk)   :: phi_L_avg, phi_R_avg ! Bulk, thickness-weighted tracer averages (left and right column)
  real, dimension(nk,2) :: phi_pp_L, phi_pp_R   ! Coefficients for the linear pseudo-reconstructions
                                                !                                                   [ nondim m^-3 ]
  real, dimension(nk)   :: h_L, h_R             ! Layer thickness (left and right)                  [m]
  real, dimension(nk)   :: khtr_u               ! Horizontal diffusivities at U-point               [m^2 s^-1]
  real                  :: hbl_L, hbl_R       ! Depth of the boundary layer (left and right)      [m]
  real                  :: F_bulk               ! Total diffusive flux across the U point           [nondim s^-1]
  real, dimension(nk)   :: F_layer              ! Diffusive flux within each layer at U-point       [nondim s^-1]
  real                  :: h_u, hblt_u          ! Thickness at the u-point                          [m]
  real                  :: khtr_avg             ! Thickness-weighted diffusivity at the u-point     [m^2 s^-1]
  real                  :: heff                 ! Harmonic mean of layer thicknesses                [m]
  real                  :: inv_heff             ! Inverse of the harmonic mean of layer thicknesses [m^[-1]
  character(len=120)    :: test_name            ! Title of the unit test
  integer               :: k_top                ! Index of cell containing top of boundary
  real                  :: zeta_top             ! Nondimension position
  integer               :: k_bot                ! Index of cell containing bottom of boundary
  real                  :: zeta_bot             ! Nondimension position

  near_boundary_unit_tests = .false.

  ! Unit tests for boundary_k_range
  test_name = 'Surface boundary spans the entire top cell'
  h_L = (/5.,5./)
  call boundary_k_range(SURFACE, nk, h_L, 5., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 1, 1., test_name, verbose)

  test_name = 'Surface boundary spans the entire column'
  h_L = (/5.,5./)
  call boundary_k_range(SURFACE, nk, h_L, 10., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 1., test_name, verbose)

  test_name = 'Bottom boundary spans the entire bottom cell'
  h_L = (/5.,5./)
  call boundary_k_range(BOTTOM, nk, h_L, 5., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 2, 0., 2, 1., test_name, verbose)

  test_name = 'Bottom boundary spans the entire column'
  h_L = (/5.,5./)
  call boundary_k_range(BOTTOM, nk, h_L, 10., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 1., test_name, verbose)

  test_name = 'Surface boundary intersects second layer'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 17.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 0.75, test_name, verbose)

  test_name = 'Surface boundary intersects first layer'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 2.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 1, 0.25, test_name, verbose)

  test_name = 'Bottom boundary intersects first layer'
  h_L = (/10.,10./)
  call boundary_k_range(BOTTOM, nk, h_L, 17.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0.75, 2, 1., test_name, verbose)

  test_name = 'Bottom boundary intersects second layer'
  h_L = (/10.,10./)
  call boundary_k_range(BOTTOM, nk, h_L, 2.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 2, 0.25, 2, 1., test_name, verbose)

  ! All cases in this section have hbl which are equal to the column thicknesses
  test_name = 'Equal hbl and same layer thicknesses (gradient from right to left)'
  hbl_L = 10; hbl_R = 10
  h_L = (/5.,5./) ; h_R = (/5.,5./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  khtr_u = (/1.,1./)
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R,&
    phi_pp_L, phi_pp_R, khtr_u, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/-5.0,-5.0/) )

  test_name = 'Equal hbl and same layer thicknesses (gradient from left to right)'
  hbl_L = 10.; hbl_R = 10.
  h_L = (/5.,5./) ; h_R = (/5.,5./)
  phi_L = (/1.,1./) ; phi_R = (/0.,0./)
  khtr_u = (/1.,1./)
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R,&
    phi_pp_L, phi_pp_R, khtr_u, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/5.0,5.0/) )

  test_name = 'Equal hbl and same layer thicknesses (no gradient)'
  hbl_L = 10; hbl_R = 10
  h_L = (/5.,5./) ; h_R = (/5.,5./)
  phi_L = (/1.,1./) ; phi_R = (/1.,1./)
  khtr_u = (/1.,1./)
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R,&
    phi_pp_L, phi_pp_R, khtr_u, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/0.0,0.0/) )

  test_name = 'Equal hbl and different layer thicknesses (gradient right to left)'
  hbl_L = 16.; hbl_R = 16.
  h_L = (/10.,6./) ; h_R = (/6.,10./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  khtr_u = (/1.,1./)
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R,&
    phi_pp_L, phi_pp_R, khtr_u, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/-8.0,-8.0/) )

  test_name = 'Equal hbl and same layer thicknesses (diagonal tracer values)'
  hbl_L = 10.; hbl_R = 10.
  h_L = (/5.,5./) ; h_R = (/5.,5./)
  phi_L = (/1.,0./) ; phi_R = (/0.,1./)
  khtr_u = (/1.,1./)
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R,&
    phi_pp_L, phi_pp_R, khtr_u, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/0.0,0.0/) )

  test_name = 'Different hbl and different column thicknesses (gradient from right to left)'
  hbl_L = 12; hbl_R = 20
  h_L = (/6.,6./) ; h_R = (/10.,10./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  khtr_u = (/1.,1./)
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R,&
    phi_pp_L, phi_pp_R, khtr_u, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/-7.5,-7.5/) )

  test_name = 'Different hbl and different layer thicknesses (gradient from right to left)'
  hbl_L = 12; hbl_R = 20
  h_L = (/6.,6./) ; h_R = (/10.,10./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  khtr_u = (/1.,1./)
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R,&
    phi_pp_L, phi_pp_R, khtr_u, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/-7.5,-7.5/) )
!
!  ! Cases where hbl < column thickness (polynomial coefficients specified for pseudo-linear reconstruction)
!  hbl_L = 2; hbl_R = 2
!  h_L = (/1.,2./) ; h_R = (/1.,2./)
!  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
!  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
!  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
!  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
!  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
!  khtr_u = (/1.,1./)
!
!  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R, phi_pp_L, phi_pp_R, khtr_u, F_layer)
!  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/-7.5,-7.5/) )


end function near_boundary_unit_tests

!> Returns true if output of near-boundary unit tests does not match correct computed values
!! and conditionally writes results to stream
logical function test_layer_fluxes(verbose, nk, test_name, F_calc, F_ans)
  logical,                    intent(in) :: verbose   !< If true, write results to stdout
  character(len=80),          intent(in) :: test_name !< Brief description of the unit test
  integer,                    intent(in) :: nk        !< Number of layers
  real, dimension(nk),        intent(in) :: F_calc    !< Fluxes of the unitless tracer from the algorithm [s^-1]
  real, dimension(nk),        intent(in) :: F_ans     !< Fluxes of the unitless tracer calculated by hand [s^-1]
  ! Local variables
  integer :: k
  integer, parameter :: stdunit = 6

  test_layer_fluxes = .false.
  do k=1,nk
    if ( F_calc(k) /= F_ans(k) ) then
      test_layer_fluxes = .true.
      write(stdunit,*) "UNIT TEST FAILED: ", test_name
      write(stdunit,10) k, F_calc(k), F_ans(k)
    elseif (verbose) then
      write(stdunit,10) k, F_calc(k), F_ans(k)
    endif
  enddo

10 format("Layer=",i3," F_calc=",f20.16," F_ans",f20.16)
end function test_layer_fluxes

!> Return true if output of unit tests for boundary_k_range does not match answers
logical function test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, k_top_ans, zeta_top_ans,&
                                       k_bot_ans, zeta_bot_ans, test_name, verbose)
  integer :: k_top               !< Index of cell containing top of boundary
  real    :: zeta_top            !< Nondimension position
  integer :: k_bot               !< Index of cell containing bottom of boundary
  real    :: zeta_bot            !< Nondimension position
  integer :: k_top_ans           !< Index of cell containing top of boundary
  real    :: zeta_top_ans        !< Nondimension position
  integer :: k_bot_ans           !< Index of cell containing bottom of boundary
  real    :: zeta_bot_ans        !< Nondimension position
  character(len=80) :: test_name !< Name of the unit test
  logical :: verbose             !< If true always print output

  integer, parameter :: stdunit = 6

  test_boundary_k_range = k_top .ne. k_top_ans
  test_boundary_k_range = test_boundary_k_range .or. (zeta_top .ne. zeta_top_ans)
  test_boundary_k_range = test_boundary_k_range .or. (k_bot .ne. k_bot_ans)
  test_boundary_k_range = test_boundary_k_range .or. (zeta_bot .ne. zeta_bot_ans)

  if (test_boundary_k_range) write(stdunit,*) "UNIT TEST FAILED: ", test_name
  if (test_boundary_k_range .or. verbose) then
    write(stdunit,20) "k_top", k_top, "k_top_ans", k_top_ans
    write(stdunit,20) "k_bot", k_bot, "k_bot_ans", k_bot_ans
    write(stdunit,30) "zeta_top", zeta_top, "zeta_top_ans", zeta_top_ans
    write(stdunit,30) "zeta_bot", zeta_bot, "zeta_bot_ans", zeta_bot_ans
  endif

  20 format(A,"=",i3,X,A,"=",i3)
  30 format(A,"=",f20.16,X,A,"=",f20.16)


end function test_boundary_k_range
end module MOM_boundary_lateral_mixing
