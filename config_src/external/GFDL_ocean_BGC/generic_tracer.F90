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
    integer, intent(in) :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes(3) !< Domain boundaries and axes
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
    real,                     intent(in) :: weight !< Unknown
    type(time_type), optional,intent(in) :: model_time !< Time
  end subroutine generic_tracer_coupler_accumulate

  !> Calls the corresponding generic_X_update_from_source routine for each package X
  subroutine generic_tracer_source(Temp,Salt,rho_dzt,dzt,hblt_depth,ilb,jlb,tau,dtts,&
       grid_dat,model_time,nbands,max_wavelength_band,sw_pen_band,opacity_band,internal_heat,&
       frunoff,grid_ht, current_wave_stress, sosga)
    real, dimension(ilb:,jlb:,:),   intent(in) :: Temp !< Potential temperature [deg C]
    real, dimension(ilb:,jlb:,:),   intent(in) :: Salt !< Salinity [psu]
    real, dimension(ilb:,jlb:,:),   intent(in) :: rho_dzt
    real, dimension(ilb:,jlb:,:),   intent(in) :: dzt !< Ocean layer thickness [m]
    real, dimension(ilb:,jlb:),     intent(in) :: hblt_depth !< Boundary layer depth
    integer,                        intent(in) :: ilb !< Lower bounds of x extent of input arrays on data domain
    integer,                        intent(in) :: jlb !< Lower bounds of y extent of input arrays on data domain
    integer,                        intent(in) :: tau !< Time step index of %field
    real,                           intent(in) :: dtts !< Unknown
    real, dimension(ilb:,jlb:),     intent(in) :: grid_dat !< Unknown
    type(time_type),                intent(in) :: model_time !< Time
    integer,                        intent(in) :: nbands !< Unknown
    real, dimension(:),             intent(in) :: max_wavelength_band
    real, dimension(:,ilb:,jlb:),   intent(in) :: sw_pen_band !< Shortwave penetration
    real, dimension(:,ilb:,jlb:,:), intent(in) :: opacity_band !< Unknown
    real, dimension(ilb:,jlb:),optional,  intent(in) :: internal_heat !< Unknown
    real, dimension(ilb:,jlb:),optional,  intent(in) :: frunoff !< Unknown
    real, dimension(ilb:,jlb:),optional,  intent(in) :: grid_ht !< Unknown
    real, dimension(ilb:,jlb:),optional , intent(in) :: current_wave_stress !< Unknown
    real,                      optional , intent(in) :: sosga ! global avg. sea surface salinity
  end subroutine generic_tracer_source

  !> Update the tracers from bottom fluxes
  subroutine generic_tracer_update_from_bottom(dt, tau, model_time)
    real,            intent(in) :: dt !< Time step increment
    integer,         intent(in) :: tau !< Time step index used for the concentration field
    type(time_type), intent(in) :: model_time !< Time
  end subroutine generic_tracer_update_from_bottom

  !> Vertically diffuse all generic tracers for GOLD ocean
  subroutine generic_tracer_vertdiff_G(h_old, ea, eb, dt, kg_m2_to_H, m_to_H, tau)
    real, dimension(:,:,:), intent(in) :: h_old !< Unknown
    real, dimension(:,:,:), intent(in) :: ea !< Unknown
    real, dimension(:,:,:), intent(in) :: eb !< Unknown
    real,                   intent(in) :: dt !< Unknown
    real,                   intent(in) :: kg_m2_to_H !< Unknown
    real,                   intent(in) :: m_to_H !< Unknown
    integer,                intent(in) :: tau !< Unknown
  end subroutine generic_tracer_vertdiff_G

  !> Set the coupler values for each generic tracer
  subroutine generic_tracer_coupler_set(IOB_struc, ST,SS,rho,ilb,jlb,tau, dzt, sosga,model_time)
    type(coupler_2d_bc_type), intent(inout) :: IOB_struc !< Ice Ocean Boundary flux structure
    integer,                     intent(in) :: ilb !< Lower bounds of x extent of input arrays on data domain
    integer,                     intent(in) :: jlb !< Lower bounds of y extent of input arrays on data domain
    integer,                     intent(in) :: tau !< Time step index of %field
    real, dimension(ilb:,jlb:),  intent(in) :: ST !< Sea surface temperature [deg C]
    real, dimension(ilb:,jlb:),  intent(in) :: SS !< Sea surface salinity [psu]
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho !< Ocean density [kg m-3]
    real, dimension(ilb:,jlb:,:), optional, intent(in) :: dzt !< Layer thickness [m]
    real,           optional, intent(in) :: sosga !< Unknown
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
