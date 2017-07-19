module coupler_indices
 
  use seq_flds_mod, only : seq_flds_x2o_fields, seq_flds_o2x_fields 
  use seq_flds_mod, only : seq_flds_i2o_per_cat, ice_ncat 
  use mct_mod

  implicit none

  private

  public coupler_indices_init

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

end module coupler_indices
