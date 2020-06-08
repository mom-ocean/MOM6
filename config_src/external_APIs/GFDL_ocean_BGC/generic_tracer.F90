module generic_tracer

  use time_manager_mod, only : time_type
  use coupler_types_mod, only : coupler_2d_bc_type

  use g_tracer_utils, only : g_tracer_type, g_diag_type

  implicit none ; private

  public generic_tracer_register
  public generic_tracer_init
  public generic_tracer_register_diag
  public generic_tracer_source
  public generic_tracer_diag
  public generic_tracer_update_from_bottom
  public generic_tracer_coupler_get
  public generic_tracer_coupler_set
  public generic_tracer_coupler_zero
  public generic_tracer_end
  public generic_tracer_get_list
  public do_generic_tracer
  public generic_tracer_vertdiff_G
  public generic_tracer_vertdiff_M
  public generic_tracer_get_diag_list
  public generic_tracer_coupler_accumulate

  logical :: do_generic_tracer = .false.

contains

  subroutine generic_tracer_register
  end subroutine generic_tracer_register

  subroutine generic_tracer_init(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes,grid_tmask,grid_kmt,init_time)
    integer,                       intent(in) :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes(3)
    type(time_type),               intent(in) :: init_time
    real, dimension(:,:,:),target, intent(in) :: grid_tmask
    integer, dimension(:,:)      , intent(in) :: grid_kmt
  end subroutine generic_tracer_init

  subroutine generic_tracer_register_diag
  end subroutine generic_tracer_register_diag

  subroutine generic_tracer_coupler_get(IOB_struc)
    type(coupler_2d_bc_type), intent(in)    :: IOB_struc
  end subroutine generic_tracer_coupler_get

  subroutine  generic_tracer_coupler_accumulate(IOB_struc, weight, model_time)
    type(coupler_2d_bc_type), intent(in)    :: IOB_struc
    real,                     intent(in)    :: weight
    type(time_type), optional,intent(in)    :: model_time
  end subroutine generic_tracer_coupler_accumulate

  subroutine generic_tracer_diag(ilb, jlb, tau, taup1, dtts, model_time, dzt, rho_dzt_tau, rho_dzt_taup1)
    integer,                        intent(in) :: ilb
    integer,                        intent(in) :: jlb
    integer,                        intent(in) :: tau
    integer,                        intent(in) :: taup1
    real,                           intent(in) :: dtts
    type(time_type),                intent(in) :: model_time
    real, dimension(ilb:,jlb:,:),   intent(in) :: dzt
    real, dimension(ilb:,jlb:,:),   intent(in) :: rho_dzt_tau
    real, dimension(ilb:,jlb:,:),   intent(in) :: rho_dzt_taup1
  end subroutine generic_tracer_diag

  subroutine generic_tracer_source(Temp,Salt,rho_dzt,dzt,hblt_depth,ilb,jlb,tau,dtts,&
       grid_dat,model_time,nbands,max_wavelength_band,sw_pen_band,opacity_band,internal_heat,&
       frunoff,grid_ht, current_wave_stress, sosga)
    real, dimension(ilb:,jlb:,:),   intent(in) :: Temp,Salt,rho_dzt,dzt
    real, dimension(ilb:,jlb:),     intent(in) :: hblt_depth
    integer,                        intent(in) :: ilb,jlb,tau
    real,                           intent(in) :: dtts
    real, dimension(ilb:,jlb:),     intent(in) :: grid_dat
    type(time_type),                intent(in) :: model_time
    integer,                        intent(in) :: nbands
    real, dimension(:),             intent(in) :: max_wavelength_band
    real, dimension(:,ilb:,jlb:),   intent(in) :: sw_pen_band
    real, dimension(:,ilb:,jlb:,:), intent(in) :: opacity_band
    real, dimension(ilb:,jlb:),optional,  intent(in) :: internal_heat
    real, dimension(ilb:,jlb:),optional,  intent(in) :: frunoff
    real, dimension(ilb:,jlb:),optional,  intent(in) :: grid_ht
    real, dimension(ilb:,jlb:),optional , intent(in) :: current_wave_stress
    real,                      optional , intent(in) :: sosga ! global avg. sea surface salinity
  end subroutine generic_tracer_source

  subroutine generic_tracer_update_from_bottom(dt, tau, model_time)
    real,    intent(in) :: dt
    integer, intent(in) ::tau
    type(time_type),                intent(in) :: model_time
  end subroutine generic_tracer_update_from_bottom

  subroutine generic_tracer_vertdiff_G(h_old, ea, eb, dt, kg_m2_to_H, m_to_H, tau)
    real, dimension(:,:,:), intent(in) :: h_old, ea, eb
    real,                   intent(in) :: dt, kg_m2_to_H, m_to_H
    integer,                intent(in) :: tau
  end subroutine generic_tracer_vertdiff_G

  subroutine generic_tracer_vertdiff_M(dh, dhw, diff_cbt, dt, Rho_0,tau)
    real, dimension(:,:,:), intent(in) :: dh, dhw, diff_cbt
    real,                   intent(in) :: dt,Rho_0
    integer,                intent(in) :: tau
  end subroutine generic_tracer_vertdiff_M

  subroutine generic_tracer_coupler_set(IOB_struc, ST,SS,rho,ilb,jlb,tau, dzt, sosga,model_time)
    type(coupler_2d_bc_type), intent(inout) :: IOB_struc
    integer, intent(in) :: ilb,jlb,tau
    real, dimension(ilb:,jlb:),  intent(in) :: ST,SS
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho
    real, dimension(ilb:,jlb:,:), optional, intent(in) :: dzt
    real,           optional, intent(in) :: sosga
    type(time_type),optional, intent(in) :: model_time
  end subroutine generic_tracer_coupler_set

  subroutine generic_tracer_coupler_zero(IOB_struc)
    type(coupler_2d_bc_type), intent(inout) :: IOB_struc
  end subroutine generic_tracer_coupler_zero

  subroutine generic_tracer_end
  end subroutine generic_tracer_end

  subroutine generic_tracer_get_list(list)
    type(g_tracer_type),    pointer    :: list
  end subroutine generic_tracer_get_list

  subroutine generic_tracer_get_diag_list(list)
    type(g_diag_type),    pointer    :: list
  end subroutine generic_tracer_get_diag_list

end module generic_tracer
