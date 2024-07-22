module cvmix_utils

!BOP
!\newpage
! !MODULE: cvmix_utils
!
! !AUTHOR:
!  Michael N. Levy, NCAR (mlevy@ucar.edu)
!
! !DESCRIPTION:
!  This module contains routines that are called by multiple modules but don't
!  specifically compute anything mixing related.
!\\
!\\

! !USES:

   use cvmix_kinds_and_types, only : cvmix_r8,                                &
                                     cvmix_strlen,                            &
                                     CVMIX_SUM_OLD_AND_NEW_VALS,              &
                                     CVMIX_MAX_OLD_AND_NEW_VALS,              &
                                     CVMIX_OVERWRITE_OLD_VAL

!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:
  public :: cvmix_update_wrap
  public :: cvmix_att_name

!EOP

contains

!BOP

! !IROUTINE: cvmix_update_wrap
! !INTERFACE:

  subroutine cvmix_update_wrap(old_vals, nlev, Mdiff_out, new_Mdiff,          &
                               Tdiff_out, new_Tdiff, Sdiff_out, new_Sdiff)

! !DESCRIPTION:
!  Update diffusivity values based on \verb|old_vals| (either overwrite, sum, or find
!  the level-by-level max)
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    integer, intent(in) :: old_vals, nlev
    real(cvmix_r8), dimension(nlev+1), optional, intent(in) :: new_Mdiff,     &
                                                               new_Tdiff,     &
                                                               new_Sdiff

! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(nlev+1), optional, intent(inout) :: Mdiff_out,  &
                                                                  Tdiff_out,  &
                                                                  Sdiff_out

!EOP
!BOC

    integer :: kw

    select case (old_vals)
      case (CVMIX_SUM_OLD_AND_NEW_VALS)
        if ((present(Mdiff_out)).and.(present(new_Mdiff))) &
          Mdiff_out = Mdiff_out + new_Mdiff
        if ((present(Tdiff_out)).and.(present(new_Tdiff))) &
          Tdiff_out = Tdiff_out + new_Tdiff
        if ((present(Sdiff_out)).and.(present(new_Sdiff))) &
          Sdiff_out = Sdiff_out + new_Sdiff
      case (CVMIX_MAX_OLD_AND_NEW_VALS)
        do kw=1,nlev+1
          if ((present(Mdiff_out)).and.(present(new_Mdiff))) &
            Mdiff_out(kw) = max(Mdiff_out(kw), new_Mdiff(kw))
          if ((present(Tdiff_out)).and.(present(new_Tdiff))) &
            Tdiff_out(kw) = max(Tdiff_out(kw), new_Tdiff(kw))
          if ((present(Sdiff_out)).and.(present(new_Sdiff))) &
            Sdiff_out(kw) = max(Sdiff_out(kw), new_Sdiff(kw))
        end do
      case (CVMIX_OVERWRITE_OLD_VAL)
        if ((present(Mdiff_out)).and.(present(new_Mdiff))) &
          Mdiff_out = new_Mdiff
        if ((present(Tdiff_out)).and.(present(new_Tdiff))) &
          Tdiff_out = new_Tdiff
        if ((present(Sdiff_out)).and.(present(new_Sdiff))) &
          Sdiff_out = new_Sdiff
      case DEFAULT
        print*, "ERROR: do not know how to handle old values!"
        stop 1
    end select

!EOC

  end subroutine cvmix_update_wrap

!BOP

! !IROUTINE: cvmix_att_name
! !INTERFACE:

  function cvmix_att_name(varname)

! !DESCRIPTION:
!  Given a variable short name, returns the precise name of the desired
!  attribute in the cvmix\_data\_type structure.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname

! !OUTPUT PARAMETERS:
    character(len=cvmix_strlen) :: cvmix_att_name

!EOP
!BOC

    select case(trim(varname))
      ! Scalars
      case ("nlev", "NumberLevels", "NumberOfLevels")
        cvmix_att_name = "nlev"
      case ("max_nlev", "MaxNumberLevels", "MaxNumberOfLevels")
        cvmix_att_name = "max_nlev"
      case ("depth", "ocn_depth", "OceanDepth", "DepthOfOcean")
        cvmix_att_name = "OceanDepth"
      case ('BoundaryLayerDepth','OBL_depth')
        cvmix_att_name = "BoundaryLayerDepth"
      case ("SSH", "surf_hgt", "SeaSurfaceHeight", "SurfaceHeight", "height")
        cvmix_att_name = "SeaSurfaceHeight"
      case ("surf_fric", "SurfaceFriction")
        cvmix_att_name = "SurfaceFriction"
      case ("surf_buoy", "SurfaceBuoyancy", "SurfaceBuoyancyForcing")
        cvmix_att_name = "SurfaceBuoyancyForcing"
      case ("lat", "latitude", "Latitude")
        cvmix_att_name = "Latitude"
      case ("lon", "longitude", "Longitude")
        cvmix_att_name = "Longitude"
      case ("coriolis", "Coriolis", "CoriolisFreq", "CoriolisFrequency")
        cvmix_att_name = "Coriolis"
      case ("kOBL_depth", "BoundaryLayerDepthIndex")
        cvmix_att_name = "kOBL_depth"
      case ("LangmuirEnhancementFactor", "EnhancementFactor",   &
              "langmuir_Efactor")
        cvmix_att_name = "LangmuirEnhancementFactor"
      case ("LangmuirNumber", "La")
        cvmix_att_name = "LangmuirNumber"
      case ("ltidal_Schmittner_socn")
        cvmix_att_name = "UseSchmittnerSouthernOceanMods"
      case ("ltidal_max")
        cvmix_att_name = "ApplyTidalMixingCap"

      ! Variables on level interfaces
      case ("zw", "zw_iface")
        cvmix_att_name = "zw_iface"
      case ("dzw", "dzw_iface")
        cvmix_att_name = "dzw"
      case ("Mdiff", "Udiff", "MomentumDiff", "MomentumDiffusivity")
        cvmix_att_name = "Mdiff_iface"
      case ("Tdiff", "TempDiff", "TemperatureDiff", "TemperatureDiffusivity")
        cvmix_att_name = "Tdiff_iface"
      case ("Sdiff", "SaltDiff", "SalinityDiff", "SalinityDiffusivity")
        cvmix_att_name = "Sdiff_iface"
      case ("Ri", "Ri_iface", "Richardson", "ShearRichardson",                &
            "RichardsonNumber", "ShearRichardsonNumber",                      &
            "ShearRichardson_iface")
        cvmix_att_name = "ShearRichardson_iface"
      case ("buoy", "buoy_iface", "N", "Nsqr", "BuoyancyFreq", "SqrBuoyancy", &
            "SqrBuoyancyFreq", "SqrBuoyancyFreq_iface")
        cvmix_att_name = "SqrBuoyancyFreq_iface"
      case ("kpp_transport", "kpp_nonlocal", "nonlocal_transport",            &
            "nonlocal", "kpp_nonlocal_iface")
        ! Note: this isn't an attribute in the data type, but put / get
        !       uses this as short hand for "both Tnonlocal and Snonlocal"
        cvmix_att_name = "kpp_nonlocal_iface"
      case ("Tnonlocal", "KPP_T_Nonlocal", "kpp_Tnonlocal", "kpp_Ttransport", &
            "kpp_Tnonlocal_iface")
        cvmix_att_name = "kpp_Tnonlocal_iface"
      case ("Snonlocal", "KPP_S_Nonlocal", "kpp_Snonlocal", "kpp_Stransport", &
            "kpp_Snonlocal_iface")
        cvmix_att_name = "kpp_Snonlocal_iface"

      ! Variables on level centers
      case ("z","zt","zt_cntr")
        cvmix_att_name = "zt_cntr"
      case ("dz", "dzt", "CellThickness")
        cvmix_att_name = "dzt"
      case ("rho", "dens", "WaterDensity", "WaterDensity_cntr")
        cvmix_att_name = "WaterDensity_cntr"
      case ("rho_lwr", "dens_lwr", "AdiabWaterDensity",                       &
            "AdiabWaterDensity_cntr")
        cvmix_att_name = "AdiabWaterDensity_cntr"
      case ("Rib", "Ri_bulk", "BulkRichardson", "BulkRichardsonNumber",       &
            "BulkRichardson_cntr")
        cvmix_att_name = "BulkRichardson_cntr"
      case ("Rrho", "strat_param")
        ! Note: this isn't an attribute in the data type, but the I/O routines
        !       use it to denote strat_param_num / strat_param_denom
        cvmix_att_name = "strat_param"
      case ("Rrho_num", "strat_param_num")
        cvmix_att_name = "strat_param_num"
      case ("Rrho_denom", "strat_param_denom")
        cvmix_att_name = "strat_param_denom"
      case ("Buoyancy","buoyancy","buoyancy_cntr")
        cvmix_att_name = "buoyancy_cntr"
      case ("U", "Vx", "Vx_cntr")
        cvmix_att_name = "Vx_cntr"
      case ("V", "Vy", "Vy_cntr")
        cvmix_att_name = "Vy_cntr"
      case ("SimmonsCoeff", "TidalCoeff")
        cvmix_att_name = "SimmonsCoeff"
      case ("SchmittnerCoeff")
        cvmix_att_name = "SchmittnerCoeff"
      case ("SchmittnerSouthernOcean")
        cvmix_att_name = "SchmittnerSouthernOcean"
      case ("exp_hab_zetar")
        cvmix_att_name = "exp_hab_zetar"
      case ("VertDep", "VertDep_iface", "vert_dep")
        cvmix_att_name = "VertDep_iface"
      case DEFAULT
        print*, "ERROR: ", trim(varname), " is not tied to an attribute of ", &
                "the cvmix_data_type structure."
        stop 1
    end select

!EOC

  end function cvmix_att_name

end module cvmix_utils
