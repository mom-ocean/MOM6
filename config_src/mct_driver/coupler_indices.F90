module coupler_indices

  ! MCT types
  use mct_mod, only : mct_aVect
  ! MCT fucntions
  use mct_mod, only : mct_avect_indexra, mct_aVect_init, mct_aVect_clean
  use seq_flds_mod, only : seq_flds_x2o_fields, seq_flds_o2x_fields
  use seq_flds_mod, only : seq_flds_i2o_per_cat, ice_ncat

  ! MOM types
  use MOM_grid,       only : ocean_grid_type
  use MOM_surface_forcing, only: ice_ocean_boundary_type
  ! MOM functions
  use MOM_domains,    only : pass_var, AGRID
  use ocean_model_mod, only : ocean_public_type

  implicit none

  private

  public coupler_indices_init
  public fill_ice_ocean_bnd
  public ocn_export

  type, public :: cpl_indices

    ! ocn -> drv

    integer :: o2x_So_t
    integer :: o2x_So_u
    integer :: o2x_So_v
    integer :: o2x_So_s
    integer :: o2x_So_dhdx
    integer :: o2x_So_dhdy
    ! QL, 150526, to wav, boundary layer depth
    integer :: o2x_So_bldepth
    integer :: o2x_Fioo_q
    integer :: o2x_Faoo_fco2_ocn
    integer :: o2x_Faoo_fdms_ocn

    ! drv -> ocn

    integer :: x2o_Si_ifrac        ! fractional ice wrt ocean
    integer :: x2o_So_duu10n       ! 10m wind speed squared           (m^2/s^2)
    integer :: x2o_Sa_pslv         ! sea-level pressure               (Pa)
    integer :: x2o_Sa_co2prog      ! bottom atm level prognostic CO2
    integer :: x2o_Sa_co2diag      ! bottom atm level diagnostic CO2
    ! QL, 150526, from wav
    integer :: x2o_Sw_lamult       ! wave model langmuir multiplier
    integer :: x2o_Sw_ustokes      ! surface Stokes drift, x-component
    integer :: x2o_Sw_vstokes      ! surface Stokes drift, y-component
    integer :: x2o_Foxx_taux       ! zonal wind stress (taux)         (W/m2   )
    integer :: x2o_Foxx_tauy       ! meridonal wind stress (tauy)     (W/m2   )
    integer :: x2o_Foxx_swnet      ! net short-wave heat flux         (W/m2   )
    integer :: x2o_Foxx_sen        ! sensible heat flux               (W/m2   )
    integer :: x2o_Foxx_lat
    integer :: x2o_Foxx_lwup       ! longwave radiation (up)          (W/m2   )
    integer :: x2o_Faxa_lwdn       ! longwave radiation (down)        (W/m2   )
    integer :: x2o_Fioi_melth      ! heat flux from snow & ice melt   (W/m2   )
    integer :: x2o_Fioi_meltw      ! snow melt flux                   (kg/m2/s)
    integer :: x2o_Fioi_bcpho      ! flux: Black Carbon hydrophobic release from sea ice component
    integer :: x2o_Fioi_bcphi      ! flux: Black Carbon hydrophilic release from sea ice component
    integer :: x2o_Fioi_flxdst     ! flux: dust release from sea ice component
    integer :: x2o_Fioi_salt       ! salt                             (kg(salt)/m2/s)
    integer :: x2o_Foxx_evap       ! evaporation flux                 (kg/m2/s)
    integer :: x2o_Faxa_prec
    integer :: x2o_Faxa_snow       ! water flux due to snow           (kg/m2/s)
    integer :: x2o_Faxa_rain       ! water flux due to rain           (kg/m2/s)
    integer :: x2o_Faxa_bcphidry   ! flux: Black   Carbon hydrophilic dry deposition
    integer :: x2o_Faxa_bcphodry   ! flux: Black   Carbon hydrophobic dry deposition
    integer :: x2o_Faxa_bcphiwet   ! flux: Black   Carbon hydrophilic wet deposition
    integer :: x2o_Faxa_ocphidry   ! flux: Organic Carbon hydrophilic dry deposition
    integer :: x2o_Faxa_ocphodry   ! flux: Organic Carbon hydrophobic dry deposition
    integer :: x2o_Faxa_ocphiwet   ! flux: Organic Carbon hydrophilic dry deposition
    integer :: x2o_Faxa_dstwet1    ! flux: Size 1 dust -- wet deposition
    integer :: x2o_Faxa_dstwet2    ! flux: Size 2 dust -- wet deposition
    integer :: x2o_Faxa_dstwet3    ! flux: Size 3 dust -- wet deposition
    integer :: x2o_Faxa_dstwet4    ! flux: Size 4 dust -- wet deposition
    integer :: x2o_Faxa_dstdry1    ! flux: Size 1 dust -- dry deposition
    integer :: x2o_Faxa_dstdry2    ! flux: Size 2 dust -- dry deposition
    integer :: x2o_Faxa_dstdry3    ! flux: Size 3 dust -- dry deposition
    integer :: x2o_Faxa_dstdry4    ! flux: Size 4 dust -- dry deposition
    integer :: x2o_Foxx_rofl       ! river runoff flux                (kg/m2/s)
    integer :: x2o_Foxx_rofi       ! ice runoff flux                  (kg/m2/s)

    ! optional per thickness category fields

    integer, dimension(:), allocatable :: x2o_frac_col      ! fraction of ocean cell, per column
    integer, dimension(:), allocatable :: x2o_fracr_col     ! fraction of ocean cell used in radiation computations, per column
    integer, dimension(:), allocatable :: x2o_qsw_fracr_col ! qsw * fracr, per column

  end type cpl_indices

contains



  subroutine coupler_indices_init(ind)

    type(cpl_indices), intent(inout) :: ind

    ! Local Variables

    type(mct_aVect) :: o2x      ! temporary
    type(mct_aVect) :: x2o      ! temporary

    integer          :: ncat  ! thickness category index
    character(len=2) :: cncat ! character version of ncat
    integer          :: ncol  ! column index
    integer          :: mcog_ncols
    integer          :: lmcog_flds_sent

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2o, rList=seq_flds_x2o_fields, lsize=1)
    call mct_aVect_init(o2x, rList=seq_flds_o2x_fields, lsize=1)

    ind%o2x_So_t          = mct_avect_indexra(o2x,'So_t')
    ind%o2x_So_u          = mct_avect_indexra(o2x,'So_u')
    ind%o2x_So_v          = mct_avect_indexra(o2x,'So_v')
    ind%o2x_So_s          = mct_avect_indexra(o2x,'So_s')
    ind%o2x_So_dhdx       = mct_avect_indexra(o2x,'So_dhdx')
    ind%o2x_So_dhdy       = mct_avect_indexra(o2x,'So_dhdy')
    ! QL, 150526, to wav, boundary layer depth
    ind%o2x_So_bldepth    = mct_avect_indexra(o2x,'So_bldepth')
    ind%o2x_Fioo_q        = mct_avect_indexra(o2x,'Fioo_q')
    ind%o2x_Faoo_fco2_ocn = mct_avect_indexra(o2x,'Faoo_fco2_ocn',perrWith='quiet')
    ind%o2x_Faoo_fdms_ocn = mct_avect_indexra(o2x,'Faoo_fdms_ocn',perrWith='quiet')
    ind%x2o_Si_ifrac      = mct_avect_indexra(x2o,'Si_ifrac')
    ind%x2o_Sa_pslv       = mct_avect_indexra(x2o,'Sa_pslv')
    ind%x2o_So_duu10n     = mct_avect_indexra(x2o,'So_duu10n')
    ! QL, 150526, from wav
    ind%x2o_Sw_lamult     = mct_avect_indexra(x2o,'Sw_lamult')
    ind%x2o_Sw_ustokes    = mct_avect_indexra(x2o,'Sw_ustokes')
    ind%x2o_Sw_vstokes    = mct_avect_indexra(x2o,'Sw_vstokes')

    ind%x2o_Foxx_tauy     = mct_avect_indexra(x2o,'Foxx_tauy')
    ind%x2o_Foxx_taux     = mct_avect_indexra(x2o,'Foxx_taux')
    ind%x2o_Foxx_swnet    = mct_avect_indexra(x2o,'Foxx_swnet')
    ind%x2o_Foxx_lat      = mct_avect_indexra(x2o,'Foxx_lat')
    ind%x2o_Foxx_sen      = mct_avect_indexra(x2o,'Foxx_sen')
    ind%x2o_Foxx_lwup     = mct_avect_indexra(x2o,'Foxx_lwup')
    ind%x2o_Faxa_lwdn     = mct_avect_indexra(x2o,'Faxa_lwdn')
    ind%x2o_Fioi_melth    = mct_avect_indexra(x2o,'Fioi_melth')
    ind%x2o_Fioi_meltw    = mct_avect_indexra(x2o,'Fioi_meltw')
    ind%x2o_Fioi_salt     = mct_avect_indexra(x2o,'Fioi_salt')
    ind%x2o_Fioi_bcpho    = mct_avect_indexra(x2o,'Fioi_bcpho')
    ind%x2o_Fioi_bcphi    = mct_avect_indexra(x2o,'Fioi_bcphi')
    ind%x2o_Fioi_flxdst   = mct_avect_indexra(x2o,'Fioi_flxdst')
    ind%x2o_Faxa_prec     = mct_avect_indexra(x2o,'Faxa_prec')
    ind%x2o_Faxa_snow     = mct_avect_indexra(x2o,'Faxa_snow')
    ind%x2o_Faxa_rain     = mct_avect_indexra(x2o,'Faxa_rain')
    ind%x2o_Foxx_evap     = mct_avect_indexra(x2o,'Foxx_evap')
    ind%x2o_Foxx_rofl     = mct_avect_indexra(x2o,'Foxx_rofl')
    ind%x2o_Foxx_rofi     = mct_avect_indexra(x2o,'Foxx_rofi')
    ind%x2o_Faxa_bcphidry = mct_avect_indexra(x2o,'Faxa_bcphidry')
    ind%x2o_Faxa_bcphodry = mct_avect_indexra(x2o,'Faxa_bcphodry')
    ind%x2o_Faxa_bcphiwet = mct_avect_indexra(x2o,'Faxa_bcphiwet')
    ind%x2o_Faxa_ocphidry = mct_avect_indexra(x2o,'Faxa_ocphidry')
    ind%x2o_Faxa_ocphodry = mct_avect_indexra(x2o,'Faxa_ocphodry')
    ind%x2o_Faxa_ocphiwet = mct_avect_indexra(x2o,'Faxa_ocphiwet')
    ind%x2o_Faxa_dstdry1  = mct_avect_indexra(x2o,'Faxa_dstdry1')
    ind%x2o_Faxa_dstdry2  = mct_avect_indexra(x2o,'Faxa_dstdry2')
    ind%x2o_Faxa_dstdry3  = mct_avect_indexra(x2o,'Faxa_dstdry3')
    ind%x2o_Faxa_dstdry4  = mct_avect_indexra(x2o,'Faxa_dstdry4')
    ind%x2o_Faxa_dstwet1  = mct_avect_indexra(x2o,'Faxa_dstwet1')
    ind%x2o_Faxa_dstwet2  = mct_avect_indexra(x2o,'Faxa_dstwet2')
    ind%x2o_Faxa_dstwet3  = mct_avect_indexra(x2o,'Faxa_dstwet3')
    ind%x2o_Faxa_dstwet4  = mct_avect_indexra(x2o,'Faxa_dstwet4')
    ind%x2o_Sa_co2prog    = mct_avect_indexra(x2o,'Sa_co2prog',perrWith='quiet')
    ind%x2o_Sa_co2diag    = mct_avect_indexra(x2o,'Sa_co2diag',perrWith='quiet')

    ! optional per thickness category fields

    ! convert cpl indices to mcog column indices
    ! this implementation only handles columns due to ice thickness categories

    lmcog_flds_sent = seq_flds_i2o_per_cat

    if (seq_flds_i2o_per_cat) then
      mcog_ncols = ice_ncat+1
      allocate(ind%x2o_frac_col(mcog_ncols))
      allocate(ind%x2o_fracr_col(mcog_ncols))
      allocate(ind%x2o_qsw_fracr_col(mcog_ncols))

      ncol = 1
      ind%x2o_frac_col(ncol)        = mct_avect_indexra(x2o,'Sf_afrac')
      ind%x2o_fracr_col(ncol)       = mct_avect_indexra(x2o,'Sf_afracr')
      ind%x2o_qsw_fracr_col(ncol)   = mct_avect_indexra(x2o,'Foxx_swnet_afracr')

      do ncat = 1, ice_ncat
        write(cncat,'(i2.2)') ncat
        ncol = ncat+1
        ind%x2o_frac_col(ncol)      = mct_avect_indexra(x2o,'Si_ifrac_'//cncat)
        ind%x2o_fracr_col(ncol)     = ind%x2o_frac_col(ncol)
        ind%x2o_qsw_fracr_col(ncol) = mct_avect_indexra(x2o,'PFioi_swpen_ifrac_'//cncat)
      enddo
    else
      mcog_ncols = 1
    endif

    call mct_aVect_clean(x2o)
    call mct_aVect_clean(o2x)


  end subroutine coupler_indices_init

  !> Maps outgoing ocean data to MCT buffer
  subroutine ocn_export(ind, ocn_public, grid, o2x)
    type(cpl_indices),       intent(inout) :: ind        !< Index
    type(ocean_public_type), intent(in)    :: ocn_public !< Ocean surface state
    type(ocean_grid_type),   intent(in)    :: grid       !< Ocean model grid
    real(kind=8),            intent(inout) :: o2x(:,:)   !< MCT outgoing bugger
    ! Local variables
    real, dimension(grid%isd:grid%ied,grid%jsd:grid%jed) :: ssh !< Local copy of sea_lev with updated halo
    integer :: i, j, n, ig, jg
    real :: slp_L, slp_R, slp_C, slope, u_min, u_max

    ! Copy from ocn_public to o2x. ocn_public uses global indexing with no halos.
    ! The mask comes from "grid" that uses the usual MOM domain that has halos
    ! and does not use global indexing.
    n = 0
    do j=grid%jsc, grid%jec
      jg = j + grid%jdg_offset
      do i=grid%isc,grid%iec
        n = n+1
        ig = i + grid%idg_offset
        o2x(ind%o2x_So_t, n) = ocn_public%t_surf(ig,jg) * grid%mask2dT(i,j)
        o2x(ind%o2x_So_s, n) = ocn_public%s_surf(ig,jg) * grid%mask2dT(i,j)
        o2x(ind%o2x_So_u, n) = ocn_public%u_surf(ig,jg) * grid%mask2dT(i,j)
        o2x(ind%o2x_So_v, n) = ocn_public%v_surf(ig,jg) * grid%mask2dT(i,j)
        ! Make a copy of ssh in order to do a halo update. We use the usual MOM domain
        ! in order to update halos. i.e. does not use global indexing.
        ssh(i,j) = ocn_public%sea_lev(ig,jg)
      end do
    end do

    ! Update halo of ssh so we can calculate gradients
    call pass_var(ssh, grid%domain)

    ! d/dx ssh
    n = 0
    do j=grid%jsc, grid%jec ; do i=grid%isc,grid%iec
      n = n+1
      ! This is a simple second-order difference
      ! o2x(ind%o2x_So_dhdx, n) = 0.5 * (ssh(i+1,j) + ssh(i-1,j)) * grid%IdxT(i,j) * grid%mask2dT(i,j)
      ! This is a PLM slope which might be less prone to the A-grid null mode
      slp_L = (ssh(i,j) - ssh(i-1,j)) * grid%mask2dCu(I-1,j)
      !if (grid%mask2dCu(I-1,j)==0.) slp_L = 0.
      slp_R = (ssh(i+1,j) - ssh(i,j)) * grid%mask2dCu(I,j)
      !if (grid%mask2dCu(I,j)==0.) slp_R = 0.
      slp_C = 0.5 * (slp_L + slp_R)
      if ( (slp_L * slp_R) > 0.0 ) then
        ! This limits the slope so that the edge values are bounded by the
        ! two cell averages spanning the edge.
        u_min = min( ssh(i-1,j), ssh(i,j), ssh(i+1,j) )
        u_max = max( ssh(i-1,j), ssh(i,j), ssh(i+1,j) )
        slope = sign( min( abs(slp_C), 2.*min( ssh(i,j) - u_min, u_max - ssh(i,j) ) ), slp_C )
      else
        ! Extrema in the mean values require a PCM reconstruction avoid generating
        ! larger extreme values.
        slope = 0.0
      end if
      o2x(ind%o2x_So_dhdx, n) = slope * grid%IdxT(i,j) * grid%mask2dT(i,j)
    end do; end do

    ! d/dy ssh
    do j=grid%jsc, grid%jec ; do i=grid%isc,grid%iec
      ! This is a simple second-order difference
    ! o2x(ind%o2x_So_dhdy, n) = 0.5 * (ssh(i,j+1) + ssh(i,j-1)) * grid%IdyT(i,j) * grid%mask2dT(i,j)
      ! This is a PLM slope which might be less prone to the A-grid null mode
      slp_L = ssh(i,j) - ssh(i,j-1)
      slp_R = ssh(i,j+1) - ssh(i,j)
slp_L=0.
slp_R=0.
      slp_C = 0.5 * (slp_L + slp_R)
      if ( (slp_L * slp_R) > 0.0 ) then
        ! This limits the slope so that the edge values are bounded by the
        ! two cell averages spanning the edge.
        u_min = min( ssh(i,j-1), ssh(i,j), ssh(i,j+1) )
        u_max = max( ssh(i,j-1), ssh(i,j), ssh(i,j+1) )
        slope = sign( min( abs(slp_C), 2.*min( ssh(i,j) - u_min, u_max - ssh(i,j) ) ), slp_C )
      else
        ! Extrema in the mean values require a PCM reconstruction avoid generating
        ! larger extreme values.
        slope = 0.0
      end if
      o2x(ind%o2x_So_dhdy, n) = slope * grid%IdyT(i,j) * grid%mask2dT(i,j)
    end do; end do

  end subroutine ocn_export


  subroutine fill_ice_ocean_bnd(ice_ocean_boundary, grid, x2o_o, ind)
    type(ice_ocean_boundary_type), intent(inout)   :: ice_ocean_boundary !< A type for the ice ocean boundary
    type(ocean_grid_type), intent(in)           :: grid
    !type(mct_aVect), intent(in)                 :: x2o_o
    real(kind=8), intent(in)                 :: x2o_o(:,:)
    type(cpl_indices), intent(inout)            :: ind

    ! local variables
    integer           :: i, j, k, ig, jg

    ! variable that are not in ice_ocean_boundary:
    ! latent (x2o_Foxx_lat)
    ! surface Stokes drift, x-comp. (x2o_Sw_ustokes)
    ! surface Stokes drift, y-comp. (x2o_Sw_vstokes)
    ! wave model langmuir multiplier (x2o_Sw_lamult)

    ! biogeochemistry
    ! Black Carbon hydrophobic release from sea ice component (x2o_Fioi_bcpho)
    ! Black Carbon hydrophilic release from sea ice component (x2o_Fioi_bcphi)
    ! dust release from sea ice component (x2o_Fioi_flxdst)
    ! Black Carbon hydrophilic dry deposition (x2o_Faxa_bcphidry)
    ! Black Carbon hydrophobic dry deposition (x2o_Faxa_bcphodry)
    ! Black Carbon hydrophobic wet deposition (x2o_Faxa_bcphiwet)
    ! Organic Carbon hydrophilic dry deposition (x2o_Faxa_ocphidry)
    ! Organic Carbon hydrophobic dry deposition (x2o_Faxa_ocphodry)
    ! Organic Carbon hydrophilic dry deposition (x2o_Faxa_ocphiwet)
    ! Sizes 1 to 4 dust - wet deposition (x2o_Faxa_dstwet?)
    ! Sizes 1 to 4 dust - dry deposition (x2o_Faxa_dstdry?)


    ! need wind_stress_multiplier?

    ! Copy from x2o to ice_ocean_boundary. ice_ocean_boundary uses global indexing with no halos.
    write(*,*) 'max. k is:', (grid%jec-grid%jsc+1) * (grid%iec-grid%isc+1)
    ! zonal wind stress (taux)
    write(*,*) 'taux', SIZE(x2o_o(ind%x2o_Foxx_taux,:))
    write(*,*) 'ice_ocean_boundary%u_flux', SIZE(ice_ocean_boundary%u_flux(:,:))
    k = 0
    do j = grid%jsc, grid%jec
      jg = j + grid%jdg_offset
      do i = grid%isc, grid%iec
        k = k + 1 ! Increment position within gindex
        ig = i + grid%idg_offset
        ! zonal wind stress (taux)
        ice_ocean_boundary%u_flux(ig,jg) = 0.0 ! x20_o(ind%x2o_Foxx_taux,k)
        ! meridional wind stress (tauy)
        ice_ocean_boundary%v_flux(ig,jg) = 0.0 ! x20_o(ind%x2o_Foxx_tauy,k)
        ! sensible heat flux
        ice_ocean_boundary%t_flux(ig,jg) = 0.0 ! x20_o(ind%x2o_Foxx_sen,k)
        ! salt flux
        ice_ocean_boundary%salt_flux(ig,jg) = 0.0 ! x20_o(ind%x2o_Fioi_salt,k)
        ! heat flux from snow & ice melt
        ice_ocean_boundary%calving_hflx(ig,jg) = 0.0 ! x20_o(ind%x2o_Fioi_melth,k)
        ! snow melt flux
        ice_ocean_boundary%fprec(ig,jg) = 0.0 ! x20_o(ind%x2o_Fioi_meltw,k)
        ! river runoff flux
        ice_ocean_boundary%runoff(ig,jg) = 0.0 ! x20_o(ind%x2o_Foxx_rofl,k)
        ! ice runoff flux
        ice_ocean_boundary%calving(ig,jg) = 0.0 ! x20_o(ind%x2o_Foxx_rofi,k)
        ! liquid precipitation (rain)
        ice_ocean_boundary%lprec(ig,jg) = 0.0 ! x20_o(ind%x2o_Faxa_rain,k)
        ! froze precipitation (snow)
        ice_ocean_boundary%fprec(ig,jg) = 0.0 ! x20_o(ind%x2o_Faxa_snow,k)
        !!!!!!! LONGWAVE NEEDS TO BE FIXED !!!!!!!
        ! longwave radiation (up)
        ice_ocean_boundary%lw_flux(ig,jg) = 0.0 ! x20_o(k,ind%x2o_Foxx_lwup)
        ! longwave radiation (down)
        ice_ocean_boundary%lw_flux(ig,jg) = 0.0 ! x20_o(k,ind%x2o_Faxa_lwdn)
        !!!!!!! SHORTWAVE NEEDS TO BE COMBINED !!!!!!!
        ! net short-wave heat flux
        ice_ocean_boundary%u_flux(ig,jg) = 0.0 ! x20_o(k,ind%x2o_Foxx_swnet)
      enddo
    enddo

    ice_ocean_boundary%wind_stagger = AGRID

  end subroutine fill_ice_ocean_bnd

end module coupler_indices
