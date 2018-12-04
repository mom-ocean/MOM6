module mom_cap_methods

  use ESMF,                only: ESMF_Clock, ESMF_ClockGet, ESMF_time, ESMF_TimeGet
  use ESMF,                only: ESMF_TimeInterval, ESMF_TimeIntervalGet
  use ESMF,                only: ESMF_State, ESMF_StateGet
  use ESMF,                only: ESMF_KIND_R8, ESMF_Field, ESMF_SUCCESS, ESMF_LogFoundError
  use ESMF,                only: ESMF_LOGERR_PASSTHRU, ESMF_FieldGet
  use ESMF,                only: ESMF_LogSetError, ESMF_RC_MEM_ALLOCATE
  use MOM_ocean_model,     only: ocean_public_type, ocean_state_type
  use MOM_surface_forcing, only: ice_ocean_boundary_type
  use MOM_grid,            only: ocean_grid_type
  use MOM_domains,         only: pass_var
  use MOM_error_handler,   only: is_root_pe
  use mpp_domains_mod,     only: mpp_get_compute_domain

  ! By default make data private
  implicit none
  private

  ! Public member functions
#ifdef CESMCOUPLED
  public :: mom_export
  public :: mom_import
#endif
  public :: mom_import_nems

  integer            :: rc,dbrc
  integer            :: import_cnt = 0
  logical, parameter :: debug=.false.

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

#ifdef CESMCOUPLED
  !> Maps outgoing ocean data to ESMF State
  !! See \ref section_mom_export for a summary of the data
  !! that is transferred from MOM6 to MCT.
  subroutine mom_export(ocean_public, grid, exportState, logunit, clock, rc)
    type(ocean_public_type) , intent(in)    :: ocean_public !< Ocean surface state
    type(ocean_grid_type)   , intent(in)    :: grid         !< Ocean model grid
    type(ESMF_State)        , intent(inout) :: exportState  !< outgoing data
    integer                 , intent(in)    :: logunit
    type(ESMF_Clock)        , intent(in)    :: clock
    integer                 , intent(inout) :: rc

    ! Local variables
    real, dimension(grid%isd:grid%ied,grid%jsd:grid%jed) :: ssh !< Local copy of sea_lev with updated halo
    integer         :: i, j, i1, j1, ig, jg, isc, iec, jsc, jec !< Grid indices
    integer         :: lbnd1, lbnd2
    real            :: slp_L, slp_R, slp_C, slope, u_min, u_max
    real            :: I_time_int  !< The inverse of coupling time interval in s-1.
    integer         :: day, secs
    type(ESMF_time) :: currTime
    real(ESMF_KIND_R8), pointer :: dataPtr_omask(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_t(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_s(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_u(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_v(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fioo_q(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_dhdx(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_dhdy(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_bldepth(:,:)
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

    lbnd1 = lbound(dataPtr_t,1)
    lbnd2 = lbound(dataPtr_t,2)

    call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

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
    do j = jsc, jec
      j1 = j + lbnd2 - jsc
      jg = j + grid%jsc - jsc
      do i = isc, iec
        i1 = i + lbnd1 - isc
        ig = i + grid%isc - isc
        dataPtr_omask(i1,j1)   = grid%mask2dT(ig,jg)
        dataPtr_t(i1,j1)       = ocean_public%t_surf(i,j) * grid%mask2dT(ig,jg) ! surface temp is in K
        dataPtr_s(i1,j1)       = ocean_public%s_surf(i,j) * grid%mask2dT(ig,jg)
        dataPtr_u(i1,j1)       = ocean_public%u_surf(i,j) * grid%mask2dT(ig,jg)
        dataPtr_v(i1,j1)       = ocean_public%v_surf(i,j) * grid%mask2dT(ig,jg)
        dataPtr_bldepth(i1,j1) = ocean_public%OBLD(i,j)   * grid%mask2dT(ig,jg)
        ! ocean melt and freeze potential (o2x_Fioo_q), W m-2
        if (ocean_public%frazil(ig,jg) > 0.0) then
           ! Frazil: change from J/m^2 to W/m^2
           dataPtr_Fioo_q(i1,j1) = ocean_public%frazil(i,j) * grid%mask2dT(ig,jg) * I_time_int
        else
           ! Melt_potential: change from J/m^2 to W/m^2
           dataPtr_Fioo_q(i1,j1) = -ocean_public%melt_potential(i,j) * grid%mask2dT(ig,jg) * I_time_int !* ncouple_per_day

           ! make sure Melt_potential is always <= 0
           if (dataPtr_Fioo_q(i1,j1) > 0.0) then
              dataPtr_Fioo_q(i1,j1) = 0.0
           endif
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
    do jg = jsc, jec
      j  = jg + grid%jsc - jsc
      j1 = jg + lbnd2 - jsc
      do ig = isc,iec
        i  = ig + grid%isc - isc
        i1 = ig + lbnd1 - isc

        ! This is a simple second-order difference
        !dataPtr_dhdx(i1,j1) = 0.5 * (ssh(i+1,j) - ssh(i-1,j)) * grid%IdxT(i,j) * grid%mask2dT(ig,jg)
        ! This is a PLM slope which might be less prone to the A-grid null mode
        slp_L = (ssh(I,j) - ssh(I-1,j)) * grid%mask2dCu(i-1,j)
        if (grid%mask2dCu(i-1,j)==0.) slp_L = 0.
        slp_R = (ssh(I+1,j) - ssh(I,j)) * grid%mask2dCu(i,j)
        if (grid%mask2dCu(i+1,j)==0.) slp_R = 0.
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
        dataPtr_dhdx(i1,j1) = slope * grid%IdxT(i,j) * grid%mask2dT(i,j)
        if (grid%mask2dT(i,j)==0.) dataPtr_dhdx(i1,j1) = 0.0
      end do
    end do

    ! d/dy ssh
    do jg = jsc, jec
      j  = jg + grid%jsc - jsc
      j1 = jg + lbnd2 - jsc
      do ig = isc,iec
        i  = ig + grid%isc - isc
        i1 = ig + lbnd1 - isc

        ! This is a simple second-order difference
        !dataPtr_dhdy(i1,j1) = 0.5 * (ssh(i,j+1) - ssh(i,j-1)) * grid%IdyT(i,j) * grid%mask2dT(ig,jg)
        ! This is a PLM slope which might be less prone to the A-grid null mode
        slp_L = ssh(i,J) - ssh(i,J-1) * grid%mask2dCv(i,j-1)
        if (grid%mask2dCv(i,j-1)==0.) slp_L = 0.
        slp_R = ssh(i,J+1) - ssh(i,J) * grid%mask2dCv(i,j)
        if (grid%mask2dCv(i,j+1)==0.) slp_R = 0.
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
        dataPtr_dhdy(i1,j1) = slope * grid%IdyT(i,j) * grid%mask2dT(i,j)
        if (grid%mask2dT(i,j)==0.) dataPtr_dhdy(i1,j1) = 0.0
      end do
    end do

    if (debug .and. is_root_pe()) then
      call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
      call ESMF_TimeGet(CurrTime, d=day, s=secs, rc=rc)

      do j = jsc, jec
        j1 = j + lbnd2 - jsc
        do i = isc, iec
          i1 = i + lbnd1 - isc
          write(logunit,F01)'export: day, secs, j, i, t_surf = ',day,secs,j,i,dataPtr_t(i1,j1)
          write(logunit,F01)'export: day, secs, j, i, s_surf = ',day,secs,j,i,dataPtr_s(i1,j1)
          write(logunit,F01)'export: day, secs, j, i, u_surf = ',day,secs,j,i,dataPtr_u(i1,j1)
          write(logunit,F01)'export: day, secs, j, i, v_surf = ',day,secs,j,i,dataPtr_v(i1,j1)
          write(logunit,F01)'export: day, secs, j, i, dhdx   = ',day,secs,j,i,dataPtr_dhdx(i1,j1)
          write(logunit,F01)'export: day, secs, j, i, dhdy   = ',day,secs,j,i,dataPtr_dhdy(i1,j1)
        end do
      end do
    end if

  end subroutine mom_export

!-----------------------------------------------------------------------

  !> This function has a few purposes: 1) it allocates and initializes the data
  !! in the fluxes structure; 2) it imports surface fluxes using data from
  !! the coupler; and 3) it can apply restoring in SST and SSS.
  !! See \ref section_ocn_import for a summary of the surface fluxes that are
  !! passed from MCT to MOM6, including fluxes that need to be included in
  !! the future.
  subroutine mom_import(ocean_public, grid, importState, ice_ocean_boundary, &
       logunit, runtype, clock, rc)

    type(ocean_public_type)       , intent(in)    :: ocean_public       !< Ocean surface state
    type(ocean_grid_type)         , intent(in)    :: grid               !< Ocean model grid
    type(ESMF_State)              , intent(inout) :: importState        !< incoming data
    type(ice_ocean_boundary_type) , intent(inout) :: ice_ocean_boundary !< Ocean boundary forcing
    type(ESMF_Clock)              , intent(in)    :: clock
    integer                       , intent(in)    :: logunit
    character(len=*)              , intent(in)    :: runtype
    integer                       , intent(inout) :: rc

    ! Local Variables
    integer                     :: i, j, i1, j1, ig, jg  ! Grid indices
    integer                     :: isc, iec, jsc, jec    ! Grid indices
    integer                     :: i0, j0, is, js, ie, je
    integer                     :: lbnd1, lbnd2
    real(ESMF_KIND_R8), pointer :: dataPtr_p(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_taux(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_tauy(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sen(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_lat(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_evap(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_osalt(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_lwdn(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_lwup(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swvdr(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swvdf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swndr(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swndf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_rofl(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_rofi(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_iosalt(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_rain(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_snow(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_lamult(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ustokes(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_vstokes(:,:)
    integer                     :: day, secs
    type(ESMF_time)             :: currTime
    logical                     :: do_import
    character(len=*), parameter :: F01  = "('(mom_import) ',a,4(i6,2x),d21.14)"
    character(len=*), parameter :: subname = '(mom_import)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call State_getFldPtr(importState,'Sa_pslv', dataPtr_p,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
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
    call State_getFldPtr(importState,"Faxa_lwdn" , dataPtr_lwdn, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"Foxx_lwup" , dataPtr_lwup, rc=rc)
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
    call State_getFldPtr(importState,"Fioi_salt" , dataPtr_iosalt, rc=rc)
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

    lbnd1 = lbound(dataPtr_p,1)
    lbnd2 = lbound(dataPtr_p,2)

    call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

    ! import_cnt is used to skip using the import state at the first count
    import_cnt = import_cnt + 1

    if ((trim(runtype) == 'initial' .and. import_cnt <= 2)) then
      ! This will skip the first time import information is given
      do_import = .false.
    else
      do_import = .true.
    end if

    if (do_import) then
       do j = jsc, jec
          j1 = j + lbnd2 - jsc
          do i = isc, iec
             i1 = i + lbnd1 - isc

             ice_ocean_boundary%p(i,j)               =  dataPtr_p(i1,j1)
             ice_ocean_boundary%u_flux(i,j)          =  dataPtr_taux(i1,j1)
             ice_ocean_boundary%v_flux(i,j)          =  dataPtr_tauy(i1,j1)
             ice_ocean_boundary%lprec(i,j)           =  dataPtr_rain(i1,j1)
             ice_ocean_boundary%fprec(i,j)           =  dataPtr_snow(i1,j1)
             ice_ocean_boundary%t_flux(i,j)          = -dataPtr_sen(i1,j1)
             ice_ocean_boundary%q_flux(i,j)          = -dataPtr_evap(i1,j1)
             ice_ocean_boundary%lw_flux(i,j)         =  dataPtr_lwup(i1,j1) + dataPtr_lwdn(i1,j1)
             ice_ocean_boundary%sw_flux_vis_dir(i,j) =  dataPtr_swvdr(i1,j1)
             ice_ocean_boundary%sw_flux_vis_dif(i,j) =  dataPtr_swvdf(i1,j1)
             ice_ocean_boundary%sw_flux_nir_dir(i,j) =  dataPtr_swndr(i1,j1)
             ice_ocean_boundary%sw_flux_nir_dif(i,j) =  dataPtr_swndf(i1,j1)
             ice_ocean_boundary%salt_flux(i,j)       =  dataPtr_iosalt(i1,j1)
             ice_ocean_boundary%runoff(i,j)          =  dataPtr_rofl(i1,j1) + dataPtr_rofi(i1,j1)
             !ice_ocean_boundary%salt_flux(i,j)      =  dataPtr_osalt(i1,j1) + ice_ocean_boundary%salt_flux(i,j)
             !ice_ocean_boundary%latent_flux(i,j)    =  dataPtr_lat(i1,j1)
          enddo
       enddo

       ! do j = jsc, jec
       !    jg = j + grid%jsc - jsc
       !    do i = isc, iec
       !       ig = i + grid%jsc - isc
       !       ice_ocean_boundary%u_flux(i,j)         = &
       !            GRID%cos_rot(ig,jg)*dataPtr_taux(i1,j1) +  GRID%sin_rot(ig,jg)*dataPtr_tauy(i1,j1)
       !       ice_ocean_boundary%v_flux(i,j)         = &
       !            GRID%cos_rot(ig,jg)*dataPtr_tauy(i1,j1) +  GRID%sin_rot(ig,jg)*dataPtr_taux(i1,j1)
       !    end do
       ! end do
    end if

    ! debug output
    if (do_import .and. debug .and. is_root_pe()) then
      call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
      call ESMF_TimeGet(CurrTime, d=day, s=secs, rc=rc)

      i0 = GRID%isc - isc
      j0 = GRID%jsc - jsc
      do j = GRID%jsc, GRID%jec
        do i = GRID%isc, GRID%iec
          write(logunit,F01)'import: day, secs, j, i, u_flux          = '&
                            ,day,secs,j,i,ice_ocean_boundary%u_flux(i-i0,j-j0)
          write(logunit,F01)'import: day, secs, j, i, v_flux          = '&
                            ,day,secs,j,i,ice_ocean_boundary%v_flux(i-i0,j-j0)
          write(logunit,F01)'import: day, secs, j, i, lprec           = '&
                            ,day,secs,j,i,ice_ocean_boundary%lprec(i-i0,j-j0)
          write(logunit,F01)'import: day, secs, j, i, lwrad           = '&
                            ,day,secs,j,i,ice_ocean_boundary%lw_flux(i-i0,j-j0)
          write(logunit,F01)'import: day, secs, j, i, q_flux          = '&
                            ,day,secs,j,i,ice_ocean_boundary%q_flux(i-i0,j-j0)
          write(logunit,F01)'import: day, secs, j, i, t_flux          = '&
                            ,day,secs,j,i,ice_ocean_boundary%t_flux(i-i0,j-j0)
          !write(logunit,F01)'import: day, secs, j, i, latent_flux     = '&
          !                  ,day,secs,j,i,ice_ocean_boundary%latent_flux(i-i0,j-j0)
          write(logunit,F01)'import: day, secs, j, i, runoff          = '&
                            ,day,secs,j,i,ice_ocean_boundary%runoff(i-i0,j-j0)
          write(logunit,F01)'import: day, secs, j, i, psurf           = '&
                            ,day,secs,j,i,ice_ocean_boundary%p(i-i0,j-j0)
          write(logunit,F01)'import: day, secs, j, i, salt_flux       = '&
                            ,day,secs,j,i,ice_ocean_boundary%salt_flux(i-i0,j-j0)
          write(logunit,F01)'import: day, secs, j, i, sw_flux_vis_dir = '&
                            ,day,secs,j,i,ice_ocean_boundary%sw_flux_vis_dir(i-i0,j-j0)
          write(logunit,F01)'import: day, secs, j, i, sw_flux_vis_dif = '&
                            ,day,secs,j,i,ice_ocean_boundary%sw_flux_vis_dif(i-i0,j-j0)
          write(logunit,F01)'import: day, secs, j, i, sw_flux_nir_dir = '&
                            ,day,secs,j,i,ice_ocean_boundary%sw_flux_nir_dir(i-i0,j-j0)
          write(logunit,F01)'import: day, secs, j, i, sw_flux_nir_dif = '&
                            ,day,secs,j,i,ice_ocean_boundary%sw_flux_nir_dir(i-i0,j-j0)
        end do
      end do
    end if

  end subroutine mom_import
#endif
  !-----------------------------------------------------------------------------

  subroutine mom_import_nems(ocean_public, grid, importState, ice_ocean_boundary, rc)

    type(ocean_public_type)       , intent(in)    :: ocean_public       !< Ocean surface state
    type(ocean_grid_type)         , intent(in)    :: grid               !< Ocean model grid
    type(ESMF_State)              , intent(inout) :: importState        !< incoming data
    type(ice_ocean_boundary_type) , intent(inout) :: ice_ocean_boundary !< Ocean boundary forcing
    integer                       , intent(inout) :: rc

    ! Local Variables
    integer                         :: i, j, i1, j1, ig, jg  ! Grid indices
    integer                         :: isc, iec, jsc, jec    ! Grid indices
    integer                         :: i0, j0, is, js, ie, je
    integer                         :: lbnd1, lbnd2
    integer                         :: ubnd1, ubnd2
    real(ESMF_KIND_R8), pointer     :: dataPtr_mask(:,:)
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

    real(ESMF_KIND_R8), allocatable :: ofld(:,:), ocz(:,:), ocm(:,:)
    real(ESMF_KIND_R8), allocatable :: mmmf(:,:), mzmf(:,:)
    integer                         :: day, secs
    type(ESMF_time)                 :: currTime
    logical                         :: do_import
    character(len=*), parameter     :: subname = '(mom_import_nems)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call State_getFldPtr(importState,"mean_zonal_moment_flx", dataPtr_mzmf, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_merid_moment_flx", dataPtr_mmmf, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_sensi_heat_flx", dataPtr_sensi, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_evap_rate" , dataPtr_evap, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_salt_rate" , dataPtr_salt, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_net_sw_ir_dif_flx" , dataPtr_swndr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,"mean_net_sw_ir_dir_flx" , dataPtr_swndf, rc=rc)
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

    lbnd1 = lbound(dataPtr_p,1)
    ubnd1 = ubound(dataPtr_p,1)
    lbnd2 = lbound(dataPtr_p,2)
    ubnd2 = ubound(dataPtr_p,2)

    call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

    allocate(mzmf(lbnd1:ubnd1,lbnd2:ubnd2))
    allocate(mmmf(lbnd1:ubnd1,lbnd2:ubnd2))
    do j  = lbnd2, ubnd2
      do i = lbnd1, ubnd1
        j1 = j - lbnd2 + jsc  ! work around local vs global indexing
        i1 = i - lbnd1 + isc
        mzmf(i,j) = grid%cos_rot(i1,j1)*dataPtr_mzmf(i,j) &
                  + grid%sin_rot(i1,j1)*dataPtr_mmmf(i,j)
        mmmf(i,j) = grid%cos_rot(i1,j1)*dataPtr_mmmf(i,j) &
                  - grid%sin_rot(i1,j1)*dataPtr_mzmf(i,j)
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
          ice_ocean_boundary%q_flux(i,j)            = -dataPtr_evap(i1,j1)
          ice_ocean_boundary%t_flux(i,j)            = -dataPtr_sensi(i1,j1)
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

  !-----------------------------------------------------------------------------

  subroutine State_GetFldPtr(ST, fldname, fldptr, rc)
    type(ESMF_State)            , intent(in)  :: ST
    character(len=*)            , intent(in)  :: fldname
    real(ESMF_KIND_R8), pointer , intent(in)  :: fldptr(:,:)
    integer, optional           , intent(out) :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(mom_cap:State_GetFldPtr)'

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
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

  end subroutine State_GetFldPtr

end module mom_cap_methods
