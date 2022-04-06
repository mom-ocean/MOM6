!> A non-functioning template of the GFDL ocean BGC
module generic_tracer

  use time_manager_mod, only : time_type
  use coupler_types_mod, only : coupler_2d_bc_type

  use g_tracer_utils, only : g_tracer_type, g_diag_type

  implicit none ; private

  public generic_tracer_register
  public generic_tracer_init
  public generic_tracer_register_diag
  public generic_tracer_source
  public generic_tracer_update_from_bottom
  public generic_tracer_coupler_get
  public generic_tracer_coupler_set
  public generic_tracer_end
  public generic_tracer_get_list
  public do_generic_tracer
  public generic_tracer_vertdiff_G
  public generic_tracer_get_diag_list
  public generic_tracer_coupler_accumulate

  !> Turn on generic tracers (note dangerous use of module data)
  logical :: do_generic_tracer = .true.

contains

  !> Unknown
  subroutine generic_tracer_register
  end subroutine generic_tracer_register

  !> Initialize generic tracers
  subroutine generic_tracer_init(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes,grid_tmask,grid_kmt,init_time)
    integer,                       intent(in) :: isc !< Computation start index in i direction
    integer,                       intent(in) :: iec !< Computation end index in i direction
    integer,                       intent(in) :: jsc !< Computation start index in j direction
    integer,                       intent(in) :: jec !< Computation end index in j direction
    integer,                       intent(in) :: isd !< Data start index in i direction
    integer,                       intent(in) :: ied !< Data end index in i direction
    integer,                       intent(in) :: jsd !< Data start index in j direction
    integer,                       intent(in) :: jed !< Data end index in j direction
    integer,                       intent(in) :: nk  !< Number of levels in k direction
    integer,                       intent(in) :: ntau !< The number of tracer time levels (always 1 for MOM6)
    integer,                       intent(in) :: axes(3) !< Domain axes?
    type(time_type),               intent(in) :: init_time !< Time
    real, dimension(:,:,:),target, intent(in) :: grid_tmask !< Mask
    integer, dimension(:,:)      , intent(in) :: grid_kmt !< Number of wet cells in column
  end subroutine generic_tracer_init

  !> Unknown
  subroutine generic_tracer_register_diag
  end subroutine generic_tracer_register_diag

  !> Get coupler values
  subroutine generic_tracer_coupler_get(IOB_struc)
    type(coupler_2d_bc_type), intent(in) :: IOB_struc !< Ice Ocean Boundary flux structure
  end subroutine generic_tracer_coupler_get

  !> Unknown
  subroutine  generic_tracer_coupler_accumulate(IOB_struc, weight, model_time)
    type(coupler_2d_bc_type), intent(in) :: IOB_struc !< Ice Ocean Boundary flux structure
    real,                     intent(in) :: weight    !< A weight for accumulating these fluxes
    type(time_type), optional,intent(in) :: model_time !< Time
  end subroutine generic_tracer_coupler_accumulate

  !> Calls the corresponding generic_X_update_from_source routine for each package X
  subroutine generic_tracer_source(Temp,Salt,rho_dzt,dzt,hblt_depth,ilb,jlb,tau,dtts,&
       grid_dat,model_time,nbands,max_wavelength_band,sw_pen_band,opacity_band,internal_heat,&
       frunoff,grid_ht, current_wave_stress, sosga)
    integer,                        intent(in) :: ilb    !< Lower bounds of x extent of input arrays on data domain
    integer,                        intent(in) :: jlb    !< Lower bounds of y extent of input arrays on data domain
    real, dimension(ilb:,jlb:,:),   intent(in) :: Temp   !< Potential temperature [deg C]
    real, dimension(ilb:,jlb:,:),   intent(in) :: Salt   !< Salinity [psu]
    real, dimension(ilb:,jlb:,:),   intent(in) :: rho_dzt !< Mass per unit area of each layer [kg m-2]
    real, dimension(ilb:,jlb:,:),   intent(in) :: dzt    !< Ocean layer thickness [m]
    real, dimension(ilb:,jlb:),     intent(in) :: hblt_depth !< Boundary layer depth [m]
    integer,                        intent(in) :: tau    !< Time step index of %field
    real,                           intent(in) :: dtts   !< The time step for this call [s]
    real, dimension(ilb:,jlb:),     intent(in) :: grid_dat !< Grid cell areas [m2]
    type(time_type),                intent(in) :: model_time !< Time
    integer,                        intent(in) :: nbands !< The number of bands of penetrating shortwave radiation
    real, dimension(:),             intent(in) :: max_wavelength_band !< The maximum wavelength in each band
                                                         !! of penetrating shortwave radiation [nm]
    real, dimension(:,ilb:,jlb:),   intent(in) :: sw_pen_band !< Penetrating shortwave radiation per band [W m-2].
                                                         !! The wavelength or angular direction band is the first index.
    real, dimension(:,ilb:,jlb:,:), intent(in) :: opacity_band !< Opacity of seawater averaged over each band [m-1].
                                                         !! The wavelength or angular direction band is the first index.
    real, dimension(ilb:,jlb:),optional,  intent(in) :: internal_heat !< Any internal or geothermal heat
                                                         !! sources that are applied to the ocean integrated
                                                         !! over this timestep [degC kg m-2]
    real, dimension(ilb:,jlb:),optional,  intent(in) :: frunoff !< Rate of iceberg calving [kg m-2 s-1]
    real, dimension(ilb:,jlb:),optional,  intent(in) :: grid_ht !< Unknown, and presently unused by MOM6
    real, dimension(ilb:,jlb:),optional , intent(in) :: current_wave_stress !< Unknown, and presently unused by MOM6
    real,                      optional , intent(in) :: sosga !< Global average sea surface salinity [ppt]
  end subroutine generic_tracer_source

  !> Update the tracers from bottom fluxes
  subroutine generic_tracer_update_from_bottom(dt, tau, model_time)
    real,            intent(in) :: dt !< Time step increment [s]
    integer,         intent(in) :: tau !< Time step index used for the concentration field
    type(time_type), intent(in) :: model_time !< Time
  end subroutine generic_tracer_update_from_bottom

  !> Vertically diffuse all generic tracers for GOLD ocean
  subroutine generic_tracer_vertdiff_G(h_old, ea, eb, dt, kg_m2_to_H, m_to_H, tau)
    real, dimension(:,:,:), intent(in) :: h_old  !< Layer thickness before entrainment [H ~> m or kg m-2]
    real, dimension(:,:,:), intent(in) :: ea     !< The amount of fluid entrained from the layer
                                                 !! above during this call [H ~> m or kg m-2]
    real, dimension(:,:,:), intent(in) :: eb     !< The amount of fluid entrained from the layer
                                                 !! below during this call [H ~> m or kg m-2]
    real,                   intent(in) :: dt     !< The amount of time covered by this call [s]
    real,                   intent(in) :: kg_m2_to_H !< A unit conversion factor from mass per unit
                                                 !! area to thickness units [H m2 kg-1 ~> m3 kg-1 or 1]
    real,                   intent(in) :: m_to_H !< A unit conversion factor from heights to
                                                 !! thickness units [H m-1 ~> 1 or kg m-3]
    integer,                intent(in) :: tau    !< The time level to work on (always 1 for MOM6)
  end subroutine generic_tracer_vertdiff_G

  !> Set the coupler values for each generic tracer
  subroutine generic_tracer_coupler_set(IOB_struc, ST,SS,rho,ilb,jlb,tau, dzt, sosga,model_time)
    type(coupler_2d_bc_type), intent(inout) :: IOB_struc !< Ice Ocean Boundary flux structure
    integer,                     intent(in) :: ilb !< Lower bounds of x extent of input arrays on data domain
    integer,                     intent(in) :: jlb !< Lower bounds of y extent of input arrays on data domain
    integer,                     intent(in) :: tau !< Time step index of %field
    real, dimension(ilb:,jlb:),  intent(in) :: ST !< Sea surface temperature [degC]
    real, dimension(ilb:,jlb:),  intent(in) :: SS !< Sea surface salinity [ppt]
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho !< Ocean density [kg m-3]
    real, dimension(ilb:,jlb:,:), optional, intent(in) :: dzt !< Layer thickness [m]
    real,           optional, intent(in) :: sosga !< Global mean sea surface salinity [ppt]
    type(time_type),optional, intent(in) :: model_time !< Time
  end subroutine generic_tracer_coupler_set

  !> End this module by calling the corresponding generic_X_end for each package X
  subroutine generic_tracer_end
  end subroutine generic_tracer_end

  !> Get a pointer to the head of the generic tracers list
  subroutine generic_tracer_get_list(list)
    type(g_tracer_type),    pointer    :: list !< Pointer to head of the linked list
  end subroutine generic_tracer_get_list

  !> Unknown
  subroutine generic_tracer_get_diag_list(list)
    type(g_diag_type),    pointer    :: list !< Pointer to head of the linked list
  end subroutine generic_tracer_get_diag_list

end module generic_tracer
