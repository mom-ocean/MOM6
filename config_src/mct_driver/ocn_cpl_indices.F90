module ocn_cpl_indices

  use mct_mod,              only: mct_avect_init, mct_avect_indexra, mct_aVect_clean, mct_aVect
  use seq_flds_mod,         only: ice_ncat, seq_flds_i2o_per_cat
  use seq_flds_mod,         only: seq_flds_x2o_fields, seq_flds_o2x_fields

  implicit none ; public

  !> Structure with indices needed for MCT attribute vectors
  type cpl_indices_type
    ! ocean to coupler
    integer :: o2x_So_t          !< Surface potential temperature (deg C)
    integer :: o2x_So_u          !< Surface zonal velocity [m s-1]
    integer :: o2x_So_v          !< Surface meridional velocity [m s-1]
    integer :: o2x_So_s          !< Surface salinity (PSU)
    integer :: o2x_So_dhdx       !< Zonal slope in the sea surface height
    integer :: o2x_So_dhdy       !< Meridional lope in the sea surface height
    integer :: o2x_So_bldepth    !< Boundary layer depth (m)
    integer :: o2x_Fioo_q        !< Ocean melt and freeze potential (W/m2)
    integer :: o2x_Faoo_fco2_ocn !< CO2 flux
    integer :: o2x_Faoo_fdms_ocn !< DMS flux

    ! coupler to ocean
    integer :: x2o_Si_ifrac      !< Fractional ice wrt ocean
    integer :: x2o_So_duu10n     !< 10m wind speed squared (m^2/s^2)
    integer :: x2o_Sa_pslv       !< Sea-level pressure (Pa)
    integer :: x2o_Sa_co2prog    !< Bottom atm level prognostic CO2
    integer :: x2o_Sa_co2diag    !< Bottom atm level diagnostic CO2
    integer :: x2o_Sw_lamult     !< Wave model langmuir multiplier
    integer :: x2o_Sw_ustokes    !< Surface Stokes drift, x-component
    integer :: x2o_Sw_vstokes    !< Surface Stokes drift, y-component
    integer :: x2o_Foxx_taux     !< Zonal wind stress (W/m2)
    integer :: x2o_Foxx_tauy     !< Meridonal wind stress (W/m2)
    integer :: x2o_Foxx_swnet    !< Net short-wave heat flux (W/m2)
    integer :: x2o_Foxx_sen      !< Sensible heat flux (W/m2)
    integer :: x2o_Foxx_lat      !< Latent heat flux  (W/m2)
    integer :: x2o_Foxx_lwup     !< Longwave radiation, up (W/m2)
    integer :: x2o_Faxa_lwdn     !< Longwave radiation, down (W/m2)
    integer :: x2o_Faxa_swvdr    !< Visible, direct shortwave (W/m2)
    integer :: x2o_Faxa_swvdf    !< Visible, diffuse shortwave (W/m2)
    integer :: x2o_Faxa_swndr    !< near-IR, direct shortwave (W/m2)
    integer :: x2o_Faxa_swndf    !< near-IR, direct shortwave (W/m2)
    integer :: x2o_Fioi_melth    !< Heat flux from snow & ice melt (W/m2)
    integer :: x2o_Fioi_meltw    !< Water flux from sea ice and snow melt (kg/m2/s)
    integer :: x2o_Fioi_bcpho    !< Black Carbon hydrophobic release from sea ice component
    integer :: x2o_Fioi_bcphi    !< Black Carbon hydrophilic release from sea ice component
    integer :: x2o_Fioi_flxdst   !< Dust release from sea ice component
    integer :: x2o_Fioi_salt     !< Salt flux    (kg(salt)/m2/s)
    integer :: x2o_Foxx_evap     !< Evaporation flux  (kg/m2/s)
    integer :: x2o_Faxa_prec     !< Total precipitation flux (kg/m2/s)
    integer :: x2o_Faxa_snow     !< Water flux due to snow (kg/m2/s)
    integer :: x2o_Faxa_rain     !< Water flux due to rain (kg/m2/s)
    integer :: x2o_Faxa_bcphidry !< Black   Carbon hydrophilic dry deposition
    integer :: x2o_Faxa_bcphodry !< Black   Carbon hydrophobic dry deposition
    integer :: x2o_Faxa_bcphiwet !< Black   Carbon hydrophilic wet deposition
    integer :: x2o_Faxa_ocphidry !< Organic Carbon hydrophilic dry deposition
    integer :: x2o_Faxa_ocphodry !< Organic Carbon hydrophobic dry deposition
    integer :: x2o_Faxa_ocphiwet !< Organic Carbon hydrophilic dry deposition
    integer :: x2o_Faxa_dstwet1  !< Size 1 dust -- wet deposition
    integer :: x2o_Faxa_dstwet2  !< Size 2 dust -- wet deposition
    integer :: x2o_Faxa_dstwet3  !< Size 3 dust -- wet deposition
    integer :: x2o_Faxa_dstwet4  !< Size 4 dust -- wet deposition
    integer :: x2o_Faxa_dstdry1  !< Size 1 dust -- dry deposition
    integer :: x2o_Faxa_dstdry2  !< Size 2 dust -- dry deposition
    integer :: x2o_Faxa_dstdry3  !< Size 3 dust -- dry deposition
    integer :: x2o_Faxa_dstdry4  !< Size 4 dust -- dry deposition
    integer :: x2o_Foxx_rofl     !< River runoff flux (kg/m2/s)
    integer :: x2o_Foxx_rofi     !< Ice runoff flux (kg/m2/s)

    ! optional per thickness category fields
    integer, dimension(:), allocatable :: x2o_frac_col      !< Fraction of ocean cell, per column
    integer, dimension(:), allocatable :: x2o_fracr_col     !< Fraction of ocean cell used  in radiation computations,
                                                            !! per column
    integer, dimension(:), allocatable :: x2o_qsw_fracr_col !< qsw * fracr, per column
  end type cpl_indices_type

  public :: cpl_indices_init

!=======================================================================
contains
!=======================================================================

  !> Determines attribute vector indices
  subroutine cpl_indices_init(ind)
    type(cpl_indices_type), intent(inout) :: ind !< Structure with coupler indices and vectors

    ! Local Variables
    type(mct_aVect)  :: o2x             !< Array with ocean to coupler data
    type(mct_aVect)  :: x2o             !< Array with coupler to ocean data
    integer          :: ncat            !< Thickness category index
    character(len=2) :: cncat           !< Character version of ncat
    integer          :: ncol            !< Column index
    integer          :: mcog_ncols      !< Number of ice thickness categories?
    integer          :: lmcog_flds_sent !< Used to convert per thickness category fields?

    ! create temporary attribute vectors
    call mct_aVect_init(x2o, rList=seq_flds_x2o_fields, lsize=1)
    call mct_aVect_init(o2x, rList=seq_flds_o2x_fields, lsize=1)

    ! ocean to coupler
    ind%o2x_So_t          = mct_avect_indexra(o2x,'So_t')
    ind%o2x_So_u          = mct_avect_indexra(o2x,'So_u')
    ind%o2x_So_v          = mct_avect_indexra(o2x,'So_v')
    ind%o2x_So_s          = mct_avect_indexra(o2x,'So_s')
    ind%o2x_So_dhdx       = mct_avect_indexra(o2x,'So_dhdx')
    ind%o2x_So_dhdy       = mct_avect_indexra(o2x,'So_dhdy')
    ind%o2x_So_bldepth    = mct_avect_indexra(o2x,'So_bldepth')
    ind%o2x_Fioo_q        = mct_avect_indexra(o2x,'Fioo_q')
    ind%o2x_Faoo_fco2_ocn = mct_avect_indexra(o2x,'Faoo_fco2_ocn',perrWith='quiet')
    ind%o2x_Faoo_fdms_ocn = mct_avect_indexra(o2x,'Faoo_fdms_ocn',perrWith='quiet')

    ! coupler to ocean
    ind%x2o_Si_ifrac      = mct_avect_indexra(x2o,'Si_ifrac')
    ind%x2o_Sa_pslv       = mct_avect_indexra(x2o,'Sa_pslv')
    ind%x2o_So_duu10n     = mct_avect_indexra(x2o,'So_duu10n')
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
    ind%x2o_Faxa_swvdr    = mct_avect_indexra(x2o,'Faxa_swvdr',perrWith='quiet')
    ind%x2o_Faxa_swvdf    = mct_avect_indexra(x2o,'Faxa_swvdf',perrWith='quiet')
    ind%x2o_Faxa_swndr    = mct_avect_indexra(x2o,'Faxa_swndr',perrWith='quiet')
    ind%x2o_Faxa_swndf    = mct_avect_indexra(x2o,'Faxa_swndf',perrWith='quiet')
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

  end subroutine cpl_indices_init

!=======================================================================

end module ocn_cpl_indices
