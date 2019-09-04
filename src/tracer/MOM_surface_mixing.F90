!> A column-wise toolbox for implementing neutral diffusion
module MOM_surface_mixing

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
use PPM_functions,             only : PPM_reconstruction, PPM_boundary_extrapolation
use regrid_edge_values,        only : edge_values_implicit_h4

implicit none ; private

#include <MOM_memory.h>
contains

!< Calculate bulk layer value of a scalar quantity as the thickness weighted average
real function bulk_average(h, hBLT, phi)
  real, dimension(:), intent(in) :: h    !< Layer thicknesses [m]
  real              , intent(in) :: hBLT !< Depth of the mixing layer [m]
  real, dimension(:), intent(in) :: phi  !< Scalar quantity
  ! Local variables
  integer :: nk ! Number of layers
  real    :: htot ! Running sum of the thicknesses (top to bottom)
  integer :: k

  ! if ( len(h) .ne. len(phi ) call MOM_error(FATAL,"surface_mixing: tracer and thicknesses of different size")
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

  harmonic_mean = (h1*h2)/(h1+h2)
end function harmonic_mean

!> Calculate the near-surface diffusive fluxes calculated from a 'bulk model'
subroutine layer_fluxes_bulk_method(nk, h_L, h_R, phi_L, phi_R, hBLT_L, hBLT_R, khtr_u, F_layer)
  integer            , intent(in   ) :: nk      !< Number of layers                        [nondim]
  real, dimension(nk), intent(in   ) :: h_L     !< Layer thickness (left)                  [m]
  real, dimension(nk), intent(in   ) :: h_R     !< Layer thickness (right)                 [m]
  real, dimension(nk), intent(in   ) :: phi_L   !< Tracer values (left)                    [ nondim m^-3 ]
  real, dimension(nk), intent(in   ) :: phi_R   !< Tracer values (right)                   [ nondim m^-3 ]
  real               , intent(in   ) :: hBLT_L  !< Depth of the boundary layer (left)      [m]
  real               , intent(in   ) :: hBLT_R  !< Depth of the boundary layer (right)     [m]
  real, dimension(nk), intent(in   ) :: khtr_u  !< Horizontal diffusivities at U-point     [m^2 s^-1]
  real, dimension(nk), intent(  out) :: F_layer !< Layerwise diffusive flux at U-point     [tracer_units s^-1]
  ! Local variables
  real                :: F_bulk               ! Total diffusive flux across the U point           [nondim s^-1]
  real, dimension(nk) :: h_means              ! Calculate the layer-wise harmonic means           [m]
  real, dimension(nk) :: h_u                  ! Thickness at the u-point                          [m]
  real                :: hblt_u               ! Boundary layer Thickness at the u-point           [m]
  real                :: khtr_avg             ! Thickness-weighted diffusivity at the u-point     [m^2 s^-1]
  real                :: heff                 ! Harmonic mean of layer thicknesses                [m]
  real                :: inv_heff             ! Inverse of the harmonic mean of layer thicknesses [m^[-1]
  real                :: phi_L_avg, phi_R_avg ! Bulk, thickness-weighted tracer averages (left and right column)
                                              !                                                   [ nondim m^-3 ]
  integer :: k
  ! Calculate bulk averages of various quantities
  phi_L_avg = bulk_average(h_L, hBLT_L, phi_L)
  phi_R_avg = bulk_average(h_R, hBLT_R, phi_R)
  do k=1,nk
    h_u(k) = 0.5 * (h_L(k) + h_R(k))
  enddo
  hblt_u = 0.5*(hBLT_L + hBLT_R)
  khtr_avg  = bulk_average(h_u, hBLT_u, khtr_u)

  ! Calculate the 'bulk' diffusive flux from the bulk averaged quantities
  heff = (hBLT_L*hBLT_R)/(hBLT_L+hBLT_R)
  F_bulk = (khtr_avg * heff) * (phi_R_avg - phi_L_avg)

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

!> Unit tests for near-surface horizontal mixing
logical function near_surface_unit_tests( verbose )
  logical,               intent(in) :: verbose !< If true, output additional information for debugging unit tests

  ! Local variables
  integer, parameter :: nk = 2                           ! Number of layers
  real, dimension(nk) :: phi_L, phi_R         ! Tracer values (left and right column)             [ nondim m^-3 ]
  real, dimension(nk) :: phi_L_avg, phi_R_avg ! Bulk, thickness-weighted tracer averages (left and right column)
                                              !                                                   [ nondim m^-3 ]
  real, dimension(nk) :: h_L, h_R             ! Layer thickness (left and right)                  [m]
  real, dimension(nk) :: khtr_u               ! Horizontal diffusivities at U-point               [m^2 s^-1]
  real                :: hBLT_L, hBLT_R       ! Depth of the boundary layer (left and right)      [m]
  real                :: F_bulk               ! Total diffusive flux across the U point           [nondim s^-1]
  real, dimension(nk) :: F_layer              ! Diffusive flux within each layer at U-point       [nondim s^-1]
  real                :: h_u, hblt_u          ! Thickness at the u-point                          [m]
  real                :: khtr_avg             ! Thickness-weighted diffusivity at the u-point     [m^2 s^-1]
  real                :: heff                 ! Harmonic mean of layer thicknesses                [m]
  real                :: inv_heff             ! Inverse of the harmonic mean of layer thicknesses [m^[-1]

  ! Equal bottom boundary layer depths and same layer thicknesses (gradient from right to left)
  hBLT_l = 10; hBLT_r = 10
  h_L = (/5,5/) ; h_R = (/5,5/)
  phi_L = (/0,0/) ; phi_R = (/1,1/)
  khtr_u = (/1,1/)

end function near_surface_unit_tests

!!> Returns true if output of find_neutral_surface_positions() does not match correct values,
!!! and conditionally writes results to stream
!logical function test_nsp(verbose, ns, KoL, KoR, pL, pR, hEff, KoL0, KoR0, pL0, pR0, hEff0, title)
!  logical,                    intent(in) :: verbose !< If true, write results to stdout
!  integer,                    intent(in) :: ns    !< Number of surfaces
!  integer, dimension(ns), intent(in) :: KoL   !< Index of first left interface above neutral surface
!  integer, dimension(ns), intent(in) :: KoR   !< Index of first right interface above neutral surface
!  real, dimension(ns),    intent(in) :: pL    !< Fractional position of neutral surface within layer KoL
!  real, dimension(ns),    intent(in) :: pR    !< Fractional position of neutral surface within layer KoR
!  real, dimension(ns-1),    intent(in) :: hEff  !< Effective thickness between two neutral surfaces [Pa]
!  integer, dimension(ns), intent(in) :: KoL0  !< Correct value for KoL
!  integer, dimension(ns), intent(in) :: KoR0  !< Correct value for KoR
!  real, dimension(ns),    intent(in) :: pL0   !< Correct value for pL
!  real, dimension(ns),    intent(in) :: pR0   !< Correct value for pR
!  real, dimension(ns-1),    intent(in) :: hEff0 !< Correct value for hEff
!  character(len=*),           intent(in) :: title !< Title for messages
!
!  ! Local variables
!  integer :: k, stdunit
!  logical :: this_row_failed
!
!  test_nsp = .false.
!  do k = 1,ns
!    test_nsp = test_nsp .or. compare_nsp_row(KoL(k), KoR(k), pL(k), pR(k), KoL0(k), KoR0(k), pL0(k), pR0(k))
!    if (k < ns) then
!      if (hEff(k) /= hEff0(k)) test_nsp = .true.
!    endif
!  enddo
!
!  if (test_nsp .or. verbose) then
!    stdunit = 6
!    if (test_nsp) stdunit = 0 ! In case of wrong results, write to error stream
!    write(stdunit,'(a)') title
!    do k = 1,ns
!      this_row_failed = compare_nsp_row(KoL(k), KoR(k), pL(k), pR(k), KoL0(k), KoR0(k), pL0(k), pR0(k))
!      if (this_row_failed) then
!        write(stdunit,10) k,KoL(k),pL(k),KoR(k),pR(k),' <-- WRONG!'
!        write(stdunit,10) k,KoL0(k),pL0(k),KoR0(k),pR0(k),' <-- should be this'
!      else
!        write(stdunit,10) k,KoL(k),pL(k),KoR(k),pR(k)
!      endif
!      if (k < ns) then
!        if (hEff(k) /= hEff0(k)) then
!          write(stdunit,'(i3,8x,"layer hEff =",2(f20.16,a))') k,hEff(k)," .neq. ",hEff0(k),' <-- WRONG!'
!        else
!          write(stdunit,'(i3,8x,"layer hEff =",f20.16)') k,hEff(k)
!        endif
!      endif
!    enddo
!  endif
!  if (test_nsp) call MOM_error(FATAL,"test_nsp failed")
!
!10 format("ks=",i3," kL=",i3," pL=",f20.16," kR=",i3," pR=",f20.16,a)
!end function test_nsp

end module MOM_surface_mixing
