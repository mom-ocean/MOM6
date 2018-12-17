module mom_cap_methods

  ! Cap import/export methods for both NEMS and CMEPS

  ! Masks, areas, center (tlat, tlon), and corner (ulat, ulon) coordinates are then added to the `ESMF_Grid`
  ! by retrieving those fields from MOM with calls to `ocean_model_data_get()`.


  use ESMF,                only: ESMF_Clock, ESMF_ClockGet, ESMF_time, ESMF_TimeGet
  use ESMF,                only: ESMF_TimeInterval, ESMF_TimeIntervalGet
  use ESMF,                only: ESMF_State, ESMF_StateGet, ESMF_Field, ESMF_FieldGet
  use ESMF,                only: ESMF_KIND_R8, ESMF_SUCCESS, ESMF_LogFoundError
  use ESMF,                only: ESMF_LOGERR_PASSTHRU, ESMF_LOGMSG_INFO, ESMF_LOGWRITE
  use ESMF,                only: ESMF_LogSetError, ESMF_RC_MEM_ALLOCATE
  use ESMF,                only: ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
  use ESMF,                only: operator(/=), operator(==)
  use MOM_ocean_model,     only: ocean_public_type, ocean_state_type, ocean_model_data_get
  use MOM_surface_forcing, only: ice_ocean_boundary_type
  use MOM_grid,            only: ocean_grid_type
  use MOM_domains,         only: pass_var
  use MOM_error_handler,   only: is_root_pe
  use mpp_domains_mod,     only: mpp_get_compute_domain

  ! By default make data private
  implicit none
  private

  ! Public member functions
  public :: mom_export_cesm
  public :: mom_import_cesm
  public :: mom_export_nems
  public :: mom_import_nems

  interface State_GetFldPtr
     module procedure State_GetFldPtr_1d
     module procedure State_GetFldPtr_2d
  end interface

  integer            :: rc,dbrc
  integer            :: import_cnt = 0

!===============================================================================
contains
!===============================================================================

  !> Maps outgoing ocean data to ESMF State
  subroutine mom_export_cesm(ocean_public, grid, exportState, logunit, clock, rc)

    ! Input/output variables
    type(ocean_public_type) , intent(in)    :: ocean_public !< Ocean surface state
    type(ocean_grid_type)   , intent(in)    :: grid         !< Ocean model grid
    type(ESMF_State)        , intent(inout) :: exportState  !< outgoing data
    integer                 , intent(in)    :: logunit
    type(ESMF_Clock)        , intent(in)    :: clock
    integer                 , intent(inout) :: rc

    ! Local variables
    real            :: ssh(grid%isd:grid%ied,grid%jsd:grid%jed) !< Local copy of sea_lev with updated halo
    integer         :: i, j, i1, j1, ig, jg, isc, iec, jsc, jec !< Grid indices
    integer         :: n
    real            :: slp_L, slp_R, slp_C, slope, u_min, u_max
    real            :: I_time_int  !< The inverse of coupling time interval in s-1.
    integer         :: day, secs
    type(ESMF_time) :: currTime
    real(ESMF_KIND_R8), pointer :: dataPtr_omask(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_t(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_s(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_u(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_v(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fioo_q(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_dhdx(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_dhdy(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_bldepth(:)
    type(ESMF_TimeInterval)     :: timeStep
    integer                     :: dt_int       !< time over which to advance the ocean (ocean_coupling_time_step), in sec
    character(len=*), parameter :: F01  = "('(mom_import) ',a,4(i6,2x),d21.14)"
    character(len=*), parameter :: subname = '(mom_export)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call State_getFldPtr(exportState,"So_omask", dataPtr_omask, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,"So_t", dataPtr_t, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,"So_s", dataPtr_s, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,"So_u", dataPtr_u, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,"So_v", dataPtr_v, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,"Fioo_q", dataPtr_fioo_q, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,"So_dhdx", dataPtr_dhdx, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,"So_dhdy", dataPtr_dhdy, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    !TODO: need to add the So_bldepth since this is needed for the wave model
    call State_getFldPtr(exportState,"So_bldepth", dataPtr_bldepth, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Use Adcroft's rule of reciprocals; it does the right thing here.
    call ESMF_ClockGet( clock, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_TimeIntervalGet( timeStep, s=dt_int, rc=rc )
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (real(dt_int) > 0.0) then
       I_time_int = 1.0 / real(dt_int)
    else
       I_time_int = 0.0
    end if

    ! Copy from ocean_public to exportstate. ocean_public uses global indexing with no halos.
    ! The mask comes from "grid" that uses the usual MOM domain that has halos
    ! and does not use global indexing.

    n = 0
    do j=grid%jsc, grid%jec
       jg = j + grid%jdg_offset
       do i=grid%isc,grid%iec
          ig = i + grid%idg_offset
          n = n+1
          dataPtr_omask(n)   = grid%mask2dT(i,j)
          dataPtr_t(n)       = ocean_public%t_surf(ig,jg) * grid%mask2dT(i,j) ! surface temp is in K
          dataPtr_s(n)       = ocean_public%s_surf(ig,jg) * grid%mask2dT(i,j)
          dataPtr_u(n)       = ocean_public%u_surf(ig,jg) * grid%mask2dT(i,j)
          dataPtr_v(n)       = ocean_public%v_surf(ig,jg) * grid%mask2dT(i,j)
          dataPtr_bldepth(n) = ocean_public%OBLD(ig,jg)   * grid%mask2dT(i,j)
          ! ocean melt and freeze potential (o2x_Fioo_q), W m-2
          if (ocean_public%frazil(ig,jg) > 0.0) then
             ! Frazil: change from J/m^2 to W/m^2
             dataPtr_Fioo_q(n) = ocean_public%frazil(ig,jg) * grid%mask2dT(i,j) * I_time_int
          else
             ! Melt_potential: change from J/m^2 to W/m^2
             dataPtr_Fioo_q(n) = -ocean_public%melt_potential(ig,jg) * grid%mask2dT(i,j) * I_time_int !* ncouple_per_day

             ! make sure Melt_potential is always <= 0
             if (dataPtr_Fioo_q(n) > 0.0) dataPtr_Fioo_q(n) = 0.0
          end if
       end do
    end do

    ! Make a copy of ssh in order to do a halo update. We use the usual MOM domain
    ! in order to update halos. i.e. does not use global indexing.
    do j=grid%jsc, grid%jec
      jg = j + grid%jdg_offset
      do i=grid%isc,grid%iec
        ig = i + grid%idg_offset
        ssh(i,j) = ocean_public%sea_lev(ig,jg)
      end do
    end do

    ! Update halo of ssh so we can calculate gradients
    call pass_var(ssh, grid%domain)

    ! d/dx ssh
    n = 0
    do j=grid%jsc, grid%jec
       do i=grid%isc,grid%iec
          n = n+1
          ! This is a simple second-order difference
          ! dataPtr_dhdx(n) = 0.5 * (ssh(i+1,j) - ssh(i-1,j)) * grid%IdxT(i,j) * grid%mask2dT(i,j)
          ! This is a PLM slope which might be less prone to the A-grid null mode
          slp_L = (ssh(I,j) - ssh(I-1,j)) * grid%mask2dCu(I-1,j)
          if (grid%mask2dCu(I-1,j)==0.) slp_L = 0.
          slp_R = (ssh(I+1,j) - ssh(I,j)) * grid%mask2dCu(I,j)
          if (grid%mask2dCu(I+1,j)==0.) slp_R = 0.
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
          endif
          dataPtr_dhdx(n) = slope * grid%IdxT(i,j) * grid%mask2dT(i,j)
          if (grid%mask2dT(i,j)==0.) dataPtr_dhdx(n) = 0.0
       enddo
    enddo

    ! d/dy ssh
    n = 0
    do j=grid%jsc, grid%jec
       do i=grid%isc,grid%iec
          n = n+1
          ! This is a simple second-order difference
          ! dataPtr_dhdy(n) = 0.5 * (ssh(i,j+1) - ssh(i,j-1)) * grid%IdyT(i,j) * grid%mask2dT(i,j)
          ! This is a PLM slope which might be less prone to the A-grid null mode
          slp_L = ssh(i,J) - ssh(i,J-1) * grid%mask2dCv(i,J-1)
          if (grid%mask2dCv(i,J-1)==0.) slp_L = 0.

          slp_R = ssh(i,J+1) - ssh(i,J) * grid%mask2dCv(i,J)
          if (grid%mask2dCv(i,J+1)==0.) slp_R = 0.

          slp_C = 0.5 * (slp_L + slp_R)
          !write(6,*)'slp_L, slp_R,i,j,slp_L*slp_R', slp_L, slp_R,i,j,slp_L*slp_R
          if ((slp_L * slp_R) > 0.0) then
             ! This limits the slope so that the edge values are bounded by the
             ! two cell averages spanning the edge.
             u_min = min( ssh(i,j-1), ssh(i,j), ssh(i,j+1) )
             u_max = max( ssh(i,j-1), ssh(i,j), ssh(i,j+1) )
             slope = sign( min( abs(slp_C), 2.*min( ssh(i,j) - u_min, u_max - ssh(i,j) ) ), slp_C )
          else
             ! Extrema in the mean values require a PCM reconstruction avoid generating
             ! larger extreme values.
             slope = 0.0
          endif
          dataPtr_dhdy(n) = slope * grid%IdyT(i,j) * grid%mask2dT(i,j)
          if (grid%mask2dT(i,j)==0.) dataPtr_dhdy(n) = 0.0
       enddo
    enddo

  end subroutine mom_export_cesm

!===============================================================================

  !> This function has a few purposes: 1) it allocates and initializes the data
  !! in the fluxes structure; 2) it imports surface fluxes using data from
  !! the coupler; and 3) it can apply restoring in SST and SSS.
  !! See \ref section_ocn_import for a summary of the surface fluxes that are
  !! passed from MCT to MOM6, including fluxes that need to be included in
  !! the future.
  subroutine mom_import_cesm(ocean_public, grid, importState, ice_ocean_boundary, &
       logunit, runtype, clock, rc)

    ! Input/output variables
    type(ocean_public_type)       , intent(in)    :: ocean_public       !< Ocean surface state
    type(ocean_grid_type)         , intent(in)    :: grid               !< Ocean model grid
    type(ESMF_State)              , intent(inout) :: importState        !< incoming data
    type(ice_ocean_boundary_type) , intent(inout) :: ice_ocean_boundary !< Ocean boundary forcing
    type(ESMF_Clock)              , intent(in)    :: clock
    integer                       , intent(in)    :: logunit
    character(len=*)              , intent(in)    :: runtype
    integer                       , intent(inout) :: rc

    ! Local Variables
    type(ESMF_StateItem_Flag)   :: itemFlag
    integer                     :: i, j, n
    integer                     :: isc, iec, jsc, jec
    integer                     :: lsize
    integer                     :: day, secs
    type(ESMF_time)             :: currTime
    logical                     :: do_import
    ! import fields that are different for cam and fv3
    logical                     :: isPresent_lwup
    logical                     :: isPresent_lwdn
    logical                     :: isPresent_lwnet
    logical                     :: isPresent_evap
    ! from atm
    real(ESMF_KIND_R8), pointer :: dataPtr_p(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_taux(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_tauy(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sen(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_lat(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_evap(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_lwdn(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_lwup(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_lwnet(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_rain(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_snow(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swvdr(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swvdf(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swndr(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swndf(:)
    ! from river
    real(ESMF_KIND_R8), pointer :: dataPtr_rofl(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_rofi(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_salt(:)
    ! from wave
    real(ESMF_KIND_R8), pointer :: dataPtr_lamult(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ustokes(:)
    real(ESMF_KIND_R8), pointer :: dataPtr_vstokes(:)
    !
    character(len=*)  , parameter :: F01  = "('(mom_import) ',a,4(i6,2x),d21.14)"
    character(len=*)  , parameter :: subname = '(mom_import)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call State_getFldPtr(importState,'Sa_pslv', dataPtr_p,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! TODO: remove these
    call State_getFldPtr(importState,"Faxa_swndr" , dataPtr_swndr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"Faxa_swndf" , dataPtr_swndf, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"Faxa_swvdr" , dataPtr_swvdr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"Faxa_swvdf" , dataPtr_swvdf, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! TODO: add these
    ! call State_getFldPtr(importState,"mean_net_sw_ir_dir_flx" , dataPtr_swndr, rc=rc)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !   line=__LINE__, &
    !   file=__FILE__)) &
    !   return  ! bail out
    ! call State_getFldPtr(importState,"mean_net_sw_ir_dif_flx" , dataPtr_swndf, rc=rc)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !   line=__LINE__, &
    !   file=__FILE__)) &
    !   return  ! bail out
    ! call State_getFldPtr(importState,"mean_net_sw_vis_dir_flx" , dataPtr_swvdr, rc=rc)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !   line=__LINE__, &
    !   file=__FILE__)) &
    !   return  ! bail out
    ! call State_getFldPtr(importState,"mean_net_sw_vis_dif_flx" , dataPtr_swvdf, rc=rc)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !   line=__LINE__, &
    !   file=__FILE__)) &
    !   return  ! bail out

    call State_getFldPtr(importState,"Foxx_taux" , dataPtr_taux, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"Foxx_tauy" , dataPtr_tauy, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"Foxx_sen"  , dataPtr_sen, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"Foxx_lat"  , dataPtr_lat, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"Foxx_evap" , dataPtr_evap, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"Foxx_rofl" , dataPtr_rofl, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"Foxx_rofi" , dataPtr_rofi, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"Fioi_salt" , dataPtr_salt, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"Faxa_rain" , dataPtr_rain, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"Faxa_snow" , dataPtr_snow, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! -------
    ! Different treatment of long wave dependent on if cam, datm or fv3
    ! -------
    ! When running with cam or datm - need Foxx_lwup and Faxa_lwdn
    ! When running with fv3 - need mean_net_lw_flx

    call ESMF_StateGet(importState, 'Foxx_lwup', itemFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
       isPresent_lwup = .true.
       call State_getFldPtr(importState,"Foxx_lwup", dataPtr_lwup, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
    else
       isPresent_lwup = .false.
    end if
    call ESMF_StateGet(importState, 'Faxa_lwdn', itemFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
       isPresent_lwdn = .true.
       call State_getFldPtr(importState, "Faxa_lwdn", dataPtr_lwdn, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
    else
       isPresent_lwdn = .false.
    end if
    call ESMF_StateGet(importState, "mean_net_lw_flx", itemFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
       isPresent_lwnet = .true.
       call State_getFldPtr(importState,"mean_net_lw_flx" , dataPtr_lwnet, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
    else
       isPresent_lwnet = .false.
    end if

    ! -------
    ! import_cnt is used to skip using the import state at the first count
    ! -------

    import_cnt = import_cnt + 1
    if ((trim(runtype) == 'initial' .and. import_cnt <= 2)) then
       do_import = .false. ! This will skip the first time import information is given
    else
       do_import = .true.
    end if

    if (do_import) then
       call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)
       n = 0
       do j = jsc,jec
          do i = isc,iec
             n = n + 1 ! Increment position within gindex
             ice_ocean_boundary%p(i,j)               =  dataPtr_p(n)     ! surface pressure
             ice_ocean_boundary%u_flux(i,j)          =  dataPtr_taux(n)  ! zonal surface stress - taux
             ice_ocean_boundary%v_flux(i,j)          =  dataPtr_tauy(n)  ! meridional surface stress - tauy
             ice_ocean_boundary%lprec(i,j)           =  dataPtr_rain(n)  ! liquid precipitation (rain)
             ice_ocean_boundary%fprec(i,j)           =  dataPtr_snow(n)  ! frozen precipitation (snow)
             ice_ocean_boundary%t_flux(i,j)          =  dataPtr_sen(n)   ! sensible heat flux (W/m2)
             ice_ocean_boundary%latent_flux(i,j)     =  dataPtr_lat(n)   ! latent heat flux (W/m^2)
             ice_ocean_boundary%q_flux(i,j)          =  dataPtr_evap(n)  ! specific humidity flux
             if (isPresent_lwup .and. isPresent_lwdn) then
                ice_ocean_boundary%lw_flux(i,j)      =  dataPtr_lwup(n) &
                                                      + dataPtr_lwdn(n)  ! longwave radiation, sum up and down (W/m2)
             else if (isPresent_lwnet) then
                ice_ocean_boundary%lw_flux(i,j)      =  dataPtr_lwnet(n) ! net longwave radiation, sum up and down (W/m2)
             end if
             ice_ocean_boundary%sw_flux_vis_dir(i,j) =  dataPtr_swvdr(n) ! visible, direct shortwave  (W/m2)
             ice_ocean_boundary%sw_flux_vis_dif(i,j) =  dataPtr_swvdf(n) ! visible, diffuse shortwave (W/m2)
             ice_ocean_boundary%sw_flux_nir_dir(i,j) =  dataPtr_swndr(n) ! near-IR, direct shortwave  (W/m2)
             ice_ocean_boundary%sw_flux_nir_dif(i,j) =  dataPtr_swndf(n) ! near-IR, diffuse shortwave (W/m2)
             ice_ocean_boundary%rofl_flux(i,j)       =  dataPtr_rofl(n)  ! ice runoff
             ice_ocean_boundary%rofi_flux(i,j)       =  dataPtr_rofi(n)  ! liquid runoff
             ice_ocean_boundary%salt_flux(i,j)       = -dataPtr_salt(n)  ! salt flux (minus sign needed here -GMM)
          enddo
       enddo
    end if

  end subroutine mom_import_cesm

!===============================================================================

  subroutine mom_export_nems(ocean_state, ocean_public, ocean_grid, dt_cpld, exportState, rc)

    ! Input/output variables
    type (ocean_state_type)  , pointer       :: ocean_state
    type (ocean_public_type) , pointer       :: ocean_public
    type (ocean_grid_type)   , pointer       :: ocean_grid
    integer                  , intent(in)    :: dt_cpld
    type(ESMF_State)         , intent(inout) :: exportState  !< outgoing data
    integer                  , intent(out)   :: rc

    ! Local variables
    integer                         :: lbnd1, lbnd2, ubnd1, ubnd2
    integer                         :: i, j, i1, j1, ig, jg !< Grid indices
    integer                         :: isc, iec, jsc, jec   !< Grid indices
    real                            :: slp_L, slp_R, slp_C, slope, u_min, u_max !JW
    real(ESMF_KIND_R8), allocatable :: ofld(:,:)
    real(ESMF_KIND_R8), allocatable :: ocz(:,:)
    real(ESMF_KIND_R8), allocatable :: ocm(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_mask(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_ocz(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_ocm(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_frazil(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_melt_potential(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_frzmlt(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_dhdx(:,:) !JW
    real(ESMF_KIND_R8), pointer     :: dataPtr_dhdy(:,:) !JW
    real(ESMF_KIND_R8), allocatable :: ssh(:,:)
    real(ESMF_KIND_R8), allocatable :: sshx(:,:)
    real(ESMF_KIND_R8), allocatable :: sshy(:,:)
    integer                         :: ijloc(2)
    character(len=240)              :: msgString
    !--------------------------------

    call State_getFldPtr(exportState,'ocn_current_zonal',dataPtr_ocz,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,'ocn_current_merid',dataPtr_ocm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    !call State_getFldPtr(exportState,'freezing_melting_potential',dataPtr_frazil,rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    ! fixfrzmlt !JW
    call State_getFldPtr(exportState,'accum_heat_frazil',dataPtr_frazil,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,'inst_melt_potential',dataPtr_melt_potential,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,'freezing_melting_potential',dataPtr_frzmlt,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,'sea_surface_slope_zonal',dataPtr_dhdx,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,'sea_surface_slope_merid',dataPtr_dhdy,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out !JW

    allocate( ssh(ocean_grid%isd:ocean_grid%ied,ocean_grid%jsd:ocean_grid%jed)) !JW
    allocate(sshx(ocean_grid%isd:ocean_grid%ied,ocean_grid%jsd:ocean_grid%jed)) !JW
    allocate(sshy(ocean_grid%isd:ocean_grid%ied,ocean_grid%jsd:ocean_grid%jed)) !JW
    ssh  = 0.0_ESMF_KIND_R8 !JW
    sshx = 0.0_ESMF_KIND_R8 !JW
    sshy = 0.0_ESMF_KIND_R8 !JW

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! note: the following code is modified from NCAR nuopc driver mom_cap_methods
    ! where is the rotation in that system?
    !
    ! Make a copy of ssh in order to do a halo update. We use the usual MOM domain
    ! in order to update halos. i.e. does not use global indexing.
    !
    ! here, isc,iec,jsc,jec are global indices on cap domain (no halos)

    do j=jsc,jec
      do i=isc,iec
        j1 = j - ocean_grid%jdg_offset
        i1 = i - ocean_grid%idg_offset
        ssh(i1,j1) = Ocean_public%sea_lev(i,j)
      end do
    end do

    ! Update halo of ssh so we can calculate gradients
    call pass_var(ssh, ocean_grid%domain)

    ! calculation of slope on native mom domains (local indexing, halos)
    ! stay inside of halos (ie 2:79,2:97)
    ! d/dx ssh
    do j = ocean_grid%jsd+1,ocean_grid%jed-1
      do i = ocean_grid%isd+1,ocean_grid%ied-1
        ! This is a simple second-order difference
        !dataPtr_dhdx(i1,j1) = 0.5 * (ssh(i+1,j) - ssh(i-1,j)) * ocean_grid%IdxT(i,j) * ocean_grid%mask2dT(ig,jg)
        ! This is a PLM slope which might be less prone to the A-grid null mode
        slp_L = (ssh(I,j) - ssh(I-1,j)) * ocean_grid%mask2dCu(i-1,j)
        if (ocean_grid%mask2dCu(i-1,j)==0.) slp_L = 0.
        slp_R = (ssh(I+1,j) - ssh(I,j)) * ocean_grid%mask2dCu(i,j)
        if (ocean_grid%mask2dCu(i+1,j)==0.) slp_R = 0.
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
        sshx(i,j) = slope * ocean_grid%IdxT(i,j) * ocean_grid%mask2dT(i,j)
        if (ocean_grid%mask2dT(i,j)==0.) sshx(i,j) = 0.0
      end do
    end do

    ! d/dy ssh
    do j = ocean_grid%jsd+1,ocean_grid%jed-1
      do i = ocean_grid%isd+1,ocean_grid%ied-1
        ! This is a simple second-order difference
        !dataPtr_dhdy(i1,j1) = 0.5 * (ssh(i,j+1) - ssh(i,j-1)) * ocean_grid%IdyT(i,j) * ocean_grid%mask2dT(ig,jg)
        ! This is a PLM slope which might be less prone to the A-grid null mode
        slp_L = ssh(i,J) - ssh(i,J-1) * ocean_grid%mask2dCv(i,j-1)
        if (ocean_grid%mask2dCv(i,j-1)==0.) slp_L = 0.
        slp_R = ssh(i,J+1) - ssh(i,J) * ocean_grid%mask2dCv(i,j)
        if (ocean_grid%mask2dCv(i,j+1)==0.) slp_R = 0.
        slp_C = 0.5 * (slp_L + slp_R)
        !write(6,*)'slp_L, slp_R,i,j,slp_L*slp_R', slp_L, slp_R,i,j,slp_L*slp_R
        if ((slp_L * slp_R) > 0.0) then
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
        sshy(i,j) = slope * ocean_grid%IdyT(i,j) * ocean_grid%mask2dT(i,j)
        if (ocean_grid%mask2dT(i,j)==0.) sshy(i,j) = 0.0
      end do
    end do
    ! rotate slopes from tripolar grid back to lat/lon grid (CCW)
    ! "grid" uses the usual MOM domain that has halos
    ! and does not use global indexing.
    ! x,y => latlon
    do j  = lbnd2, ubnd2
      do i = lbnd1, ubnd1
        j1 = j + ocean_grid%jsc - lbnd2
        i1 = i + ocean_grid%isc - lbnd1
        dataPtr_dhdx(i,j) = ocean_grid%cos_rot(i1,j1)*sshx(i1,j1) &
                          + ocean_grid%sin_rot(i1,j1)*sshy(i1,j1)
        dataPtr_dhdy(i,j) = ocean_grid%cos_rot(i1,j1)*sshy(i1,j1) &
                          - ocean_grid%sin_rot(i1,j1)*sshx(i1,j1)
      enddo
    enddo
    deallocate(ssh); deallocate(sshx); deallocate(sshy)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    dataPtr_frazil = dataPtr_frazil/dt_cpld !convert from J/m^2 to W/m^2 for CICE coupling

    dataPtr_melt_potential = -dataPtr_melt_potential/dt_cpld !convert from J/m^2 to W/m^2 for CICE coupling
                                                             !melt_potential, defined positive for T>Tfreeze
                                                             !so change sign
    !testing
    ijloc = maxloc(dataPtr_frazil)
    if((sum(ijloc) .gt. 2) .and. (dataPtr_frazil(ijloc(1),ijloc(2)) .gt. 0.0))then
       i1 = ijloc(1) - lbnd1 + isc
       j1 = ijloc(2) - lbnd2 + jsc  ! work around local vs global indexing

       write (msgString,*)' MOM6 dataPtr_frazil at maxloc ',i1,j1,&
            real(dataPtr_frazil(ijloc(1),ijloc(2)),4)
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

       write (msgString,*)' MOM6 dataPtr_melt_potential at maxloc ',i1,j1,&
            real(dataPtr_melt_potential(ijloc(1),ijloc(2)),4)
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    endif
    !testing

    dataPtr_melt_potential = min(dataPtr_melt_potential,0.0)

    do j  = lbnd2, ubnd2
      do i = lbnd1, ubnd1
       if(dataPtr_frazil(i,j) .eq. 0.0)then
        dataPtr_frzmlt(i,j) = dataPtr_melt_potential(i,j)
       else
        dataPtr_frzmlt(i,j) = dataPtr_frazil(i,j)
       endif
      enddo
    enddo
    dataPtr_frzmlt = max(-1000.0,min(1000.0,dataPtr_frzmlt))

    ! rotate ocn current from tripolar grid back to lat/lon grid (CCW)
    ! "grid" uses the usual MOM domain that has halos and does not use global indexing.
    ! x,y => latlon

    allocate(ofld(isc:iec,jsc:jec))
    call ocean_model_data_get(ocean_state, ocean_public, 'mask', ofld, isc, jsc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
      j1 = j - lbnd2 + jsc
      i1 = i - lbnd1 + isc
      dataPtr_mask(i,j) = nint(ofld(i1,j1))
    enddo
    enddo
    deallocate(ofld)

    allocate(ocz(lbnd1:ubnd1,lbnd2:ubnd2))
    allocate(ocm(lbnd1:ubnd1,lbnd2:ubnd2))
    ocz = dataPtr_ocz
    ocm = dataPtr_ocm
    do j  = lbnd2, ubnd2
      do i = lbnd1, ubnd1
        j1 = j + ocean_grid%jsc - lbnd2
        i1 = i + ocean_grid%isc - lbnd1
        dataPtr_ocz(i,j) = ocean_grid%cos_rot(i1,j1)*ocz(i,j) &
                         + ocean_grid%sin_rot(i1,j1)*ocm(i,j)
        dataPtr_ocm(i,j) = ocean_grid%cos_rot(i1,j1)*ocm(i,j) &
                         - ocean_grid%sin_rot(i1,j1)*ocz(i,j)
        ! multiply by mask to zero out non-ocean points
        dataPtr_ocz(i,j) = dataPtr_ocz(i,j) * dataPtr_mask(i,j)
        dataPtr_ocm(i,j) = dataPtr_ocm(i,j) * dataPtr_mask(i,j)
      enddo
    enddo
    deallocate(ocz, ocm)

  end subroutine mom_export_nems

!===============================================================================

  subroutine mom_import_nems(ocean_public, ocean_grid, importState, ice_ocean_boundary, rc)

    ! Input/output variables
    type(ocean_public_type)       , intent(in)    :: ocean_public       !< Ocean surface state
    type(ocean_grid_type)         , intent(in)    :: ocean_grid         !< Ocean model grid
    type(ESMF_State)              , intent(inout) :: importState        !< incoming data
    type(ice_ocean_boundary_type) , intent(inout) :: ice_ocean_boundary !< Ocean boundary forcing
    integer                       , intent(inout) :: rc

    ! Local Variables
    integer                         :: lbnd1,ubnd1,lbnd2,ubnd2
    integer                         :: i, j, i1, j1, ig, jg   ! Grid indices
    integer                         :: isc, iec, jsc, jec     ! Grid indices
    integer                         :: i0, j0, is, js, ie, je ! Grid indices
    real(ESMF_KIND_R8), pointer     :: dataPtr_p(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_mmmf(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_mzmf(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_sensi(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_evap(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_salt(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_lwflux(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_swvdr(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_swvdf(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_swndr(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_swndf(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_runoff(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_rain(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_snow(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_calving(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_runoff_hflx(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_calving_hflx(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtr_mi(:,:)

    real(ESMF_KIND_R8), allocatable :: mmmf(:,:), mzmf(:,:)
    integer                         :: day, secs
    type(ESMF_time)                 :: currTime
    logical                         :: do_import
    character(len=*), parameter     :: subname = '(mom_import_nems)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call State_getFldPtr(importState,'mean_zonal_moment_flx',dataPtr_mzmf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,'mean_merid_moment_flx',dataPtr_mmmf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,'mean_evap_rate',dataPtr_evap,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,'mean_sensi_heat_flx',dataPtr_sensi,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_salt_rate" , dataPtr_salt, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_net_sw_ir_dir_flx" , dataPtr_swndr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_net_sw_ir_dif_flx" , dataPtr_swndf, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_net_sw_vis_dir_flx" , dataPtr_swvdr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_net_sw_vis_dif_flx" , dataPtr_swvdf, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_prec_rate" , dataPtr_rain, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_fprec_rate" , dataPtr_snow, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_runoff_rate" , dataPtr_runoff, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_calving_rate" , dataPtr_calving, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_runoff_heat_flux" , dataPtr_runoff_hflx, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_calving_heat_flux" , dataPtr_calving_hflx, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,'inst_pres_height_surface', dataPtr_p,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mass_of_overlying_ice" , dataPtr_mi, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

    lbnd1 = lbound(dataPtr_p,1)
    ubnd1 = ubound(dataPtr_p,1)
    lbnd2 = lbound(dataPtr_p,2)
    ubnd2 = ubound(dataPtr_p,2)
    print *, 'lbnd1,ubnd1,lbnd2,ubnd2', lbnd1, ubnd1, lbnd2, ubnd2

    allocate(mzmf(lbnd1:ubnd1,lbnd2:ubnd2))
    allocate(mmmf(lbnd1:ubnd1,lbnd2:ubnd2))
    do j  = lbnd2, ubnd2
      do i = lbnd1, ubnd1
        j1 = j + ocean_grid%jsc - lbnd2
        i1 = i + ocean_grid%isc - lbnd1
        mzmf(i,j) = ocean_grid%cos_rot(i1,j1)*dataPtr_mzmf(i,j) &
                  - ocean_grid%sin_rot(i1,j1)*dataPtr_mmmf(i,j)
        mmmf(i,j) = ocean_grid%cos_rot(i1,j1)*dataPtr_mmmf(i,j) &
                  + ocean_grid%sin_rot(i1,j1)*dataPtr_mzmf(i,j)
      enddo
    enddo
    dataPtr_mzmf = mzmf
    dataPtr_mmmf = mmmf
    deallocate(mzmf, mmmf)

    do j = jsc, jec
       j1 = j + lbnd2 - jsc
       do i = isc, iec
          i1 = i + lbnd1 - isc

          ice_ocean_boundary%u_flux(i,j)            =  dataPtr_mzmf(i1,j1)
          ice_ocean_boundary%v_flux(i,j)            =  dataPtr_mmmf(i1,j1)
          ice_ocean_boundary%q_flux(i,j)            =  dataPtr_evap(i1,j1)
          ice_ocean_boundary%t_flux(i,j)            =  dataPtr_sensi(i1,j1)
          ice_ocean_boundary%salt_flux(i,j)         =  dataPtr_salt(i1,j1)
          ice_ocean_boundary%lw_flux(i,j)           =  dataPtr_lwflux(i1,j1)
          ice_ocean_boundary%sw_flux_vis_dir(i,j)   =  dataPtr_swvdr(i1,j1)
          ice_ocean_boundary%sw_flux_vis_dif(i,j)   =  dataPtr_swvdf(i1,j1)
          ice_ocean_boundary%sw_flux_nir_dir(i,j)   =  dataPtr_swndr(i1,j1)
          ice_ocean_boundary%sw_flux_nir_dif(i,j)   =  dataPtr_swndf(i1,j1)
          ice_ocean_boundary%lprec(i,j)             =  dataPtr_rain(i1,j1)
          ice_ocean_boundary%fprec(i,j)             =  dataPtr_snow(i1,j1)
          ice_ocean_boundary%runoff(i,j)            =  dataPtr_runoff(i1,j1)
          ice_ocean_boundary%calving(i,j)           =  dataPtr_calving(i1,j1)
          ice_ocean_boundary%runoff_hflx(i,j)       =  dataPtr_runoff_hflx(i1,j1)
          ice_ocean_boundary%calving_hflx(i,j)      =  dataPtr_calving_hflx(i1,j1)
          ice_ocean_boundary%p(i,j)                 =  dataPtr_p(i1,j1)
          ice_ocean_boundary%mi(i,j)                =  dataPtr_mi(i1,j1)
       enddo
    enddo

  end subroutine mom_import_nems

!===============================================================================

  subroutine State_GetFldPtr_1d(State, fldname, fldptr, rc)
    type(ESMF_State)            , intent(in)  :: State
    character(len=*)            , intent(in)  :: fldname
    real(ESMF_KIND_R8), pointer , intent(in)  :: fldptr(:)
    integer, optional           , intent(out) :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(mom_cap:State_GetFldPtr)'

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (present(rc)) rc = lrc

  end subroutine State_GetFldPtr_1d

!===============================================================================

  subroutine State_GetFldPtr_2d(State, fldname, fldptr, rc)
    type(ESMF_State)            , intent(in)  :: State
    character(len=*)            , intent(in)  :: fldname
    real(ESMF_KIND_R8), pointer , intent(in)  :: fldptr(:,:)
    integer, optional           , intent(out) :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(mom_cap:State_GetFldPtr)'

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (present(rc)) rc = lrc

  end subroutine State_GetFldPtr_2d

end module mom_cap_methods
