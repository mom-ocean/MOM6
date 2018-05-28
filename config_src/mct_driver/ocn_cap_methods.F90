module ocn_cap_methods

  use MOM_ocean_model,     only: ocean_public_type, ocean_state_type
  use MOM_surface_forcing, only: ice_ocean_boundary_type
  use MOM_grid,            only: ocean_grid_type
  use MOM_domains,         only: pass_var
  use mpp_domains_mod,     only: mpp_get_compute_domain
  use ocn_cpl_indices,     only: cpl_indices_type

  implicit none
  private

  public :: ocn_import
  public :: ocn_export

!=======================================================================
contains
!=======================================================================

  subroutine ocn_import(x2o, ind, grid, ice_ocean_boundary, c1, c2, c3, c4)
    real(kind=8)                  , intent(in)    :: x2o(:,:)           !< incoming data
    type(cpl_indices_type)        , intent(in)    :: ind                !< Structure with MCT attribute vectors and indices
    type(ocean_grid_type)         , intent(in)    :: grid               !< Ocean model grid
    type(ice_ocean_boundary_type) , intent(inout) :: ice_ocean_boundary !< Ocean boundary forcing
    real(kind=8), optional        , intent(in)    :: c1, c2, c3, c4     !< Coeffs. used in the shortwave decomposition

    ! Local variables
    integer :: i, j, i1, j1, ig, jg, isc, iec, jsc, jec  !< Grid indices
    integer :: k !< temporary
    !-----------------------------------------------------------------------

    isc = GRID%isc; iec = GRID%iec ; jsc = GRID%jsc; jec = GRID%jec

    k = 0
    do j = jsc, jec
       jg = j + grid%jsc - jsc
       do i = isc, iec
          ig = i + grid%jsc - isc
          k = k + 1 ! Increment position within gindex

          ! liquid precipitation (rain)
          ice_ocean_boundary%lprec(i,j) = x2o(ind%x2o_Faxa_rain,k) * GRID%mask2dT(ig,jg)

          ! frozen precipitation (snow)
          ice_ocean_boundary%fprec(i,j) = x2o(ind%x2o_Faxa_snow,k) * GRID%mask2dT(ig,jg)

          ! longwave radiation, sum up and down (W/m2)
          ice_ocean_boundary%lw_flux(i,j) = (x2o(ind%x2o_Faxa_lwdn,k) + x2o(ind%x2o_Foxx_lwup,k)) * GRID%mask2dT(i,j)

          ! sensible heat flux (W/m2)
          ice_ocean_boundary%t_flux(i,j) = x2o(ind%x2o_Foxx_sen,k) * GRID%mask2dT(i,j)

          ! latent heat flux (W/m^2)
          ice_ocean_boundary%latent_flux(i,j) = x2o(ind%x2o_Foxx_lat,k) * GRID%mask2dT(i,j)

          ! 1) visible, direct shortwave  (W/m2)
          ! 2) visible, diffuse shortwave (W/m2)
          ! 3) near-IR, direct shortwave  (W/m2)
          ! 4) near-IR, diffuse shortwave (W/m2)
          if (present(c1) .and. present(c2) .and. present(c3) .and. present(c4)) then
             ! Use runtime coefficients to decompose net short-wave heat flux into 4 components
             ice_ocean_boundary%sw_flux_vis_dir(i,j) = x2o(ind%x2o_Foxx_swnet,k) * c1 * GRID%mask2dT(ig,jg)
             ice_ocean_boundary%sw_flux_vis_dif(i,j) = x2o(ind%x2o_Foxx_swnet,k) * c2 * GRID%mask2dT(ig,jg)
             ice_ocean_boundary%sw_flux_nir_dir(i,j) = x2o(ind%x2o_Foxx_swnet,k) * c3 * GRID%mask2dT(ig,jg)
             ice_ocean_boundary%sw_flux_nir_dif(i,j) = x2o(ind%x2o_Foxx_swnet,k) * c4 * GRID%mask2dT(ig,jg)
          else
             ice_ocean_boundary%sw_flux_vis_dir(i,j) = x2o(ind%x2o_Faxa_swvdr,k) * GRID%mask2dT(ig,jg)
             ice_ocean_boundary%sw_flux_vis_dif(i,j) = x2o(ind%x2o_Faxa_swvdf,k) * GRID%mask2dT(ig,jg)
             ice_ocean_boundary%sw_flux_nir_dir(i,j) = x2o(ind%x2o_Faxa_swndr,k) * GRID%mask2dT(ig,jg)
             ice_ocean_boundary%sw_flux_nir_dif(i,j) = x2o(ind%x2o_Faxa_swndf,k) * GRID%mask2dT(ig,jg)
          end if
       end do
    end do

  end subroutine ocn_import

!=======================================================================

  !> Maps outgoing ocean data to MCT attribute vector real array
  subroutine ocn_export(ocean_public, grid, o2x, rc)
    type(ocean_public_type) , intent(in)    :: ocean_public !< Ocean surface state
    type(ocean_grid_type)   , intent(in)    :: grid         !< Ocean model grid
    real(kind=8)            , intent(out)   :: o2x(:,:)     !< outgoing data
    integer                 , intent(inout) :: rc 

    ! Local variables
    real, dimension(grid%isd:grid%ied,grid%jsd:grid%jed) :: ssh !< Local copy of sea_lev with updated halo
    integer :: i, j, i1, j1, ig, jg, isc, iec, jsc, jec         !< Grid indices
    integer :: lbnd1, lbnd2, ubnd1, ubnd2
    real :: slp_L, slp_R, slp_C, slope, u_min, u_max
    !-----------------------------------------------------------------------
    ! Nothing for now
  end subroutine ocn_export

end module ocn_cap_methods
