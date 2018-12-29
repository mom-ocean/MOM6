module mom_cap_methods

  ! Cap import/export methods for both NEMS and CMEPS

  ! Masks, areas, center (tlat, tlon), and corner (ulat, ulon) coordinates are then added to the `ESMF_Grid`
  ! by retrieving those fields from MOM with calls to `ocean_model_data_get()`.


  use NUOPC,               only: NUOPC_Advertise, NUOPC_Realize, NUOPC_IsConnected
  use NUOPC_Model,         only: NUOPC_ModelGet
  use ESMF,                only: ESMF_Clock, ESMF_ClockGet, ESMF_time, ESMF_TimeGet
  use ESMF,                only: ESMF_TimeInterval, ESMF_TimeIntervalGet
  use ESMF,                only: ESMF_State, ESMF_StateGet, ESMF_StateRemove
  use ESMF,                only: ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
  use ESMF,                only: ESMF_GridComp, ESMF_Mesh, ESMF_Grid, ESMF_GridCreate
  use ESMF,                only: ESMF_DistGrid, ESMF_DistGridCreate
  use ESMF,                only: ESMF_KIND_R8, ESMF_SUCCESS, ESMF_LogFoundError
  use ESMF,                only: ESMF_LOGERR_PASSTHRU, ESMF_LOGMSG_INFO, ESMF_LOGWRITE
  use ESMF,                only: ESMF_LogSetError, ESMF_RC_MEM_ALLOCATE
  use ESMF,                only: ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
  use ESMF,                only: ESMF_GEOMTYPE_FLAG, ESMF_GEOMTYPE_GRID, ESMF_GEOMTYPE_MESH
  use ESMF,                only: ESMF_RC_VAL_OUTOFRANGE, ESMF_INDEX_DELOCAL, ESMF_MESHLOC_ELEMENT
  use ESMF,                only: ESMF_TYPEKIND_R8
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
  public :: mom_import
  public :: mom_export_cesm
  public :: mom_export_nems

  private :: state_getimport

  interface State_GetFldPtr
     module procedure State_GetFldPtr_1d
     module procedure State_GetFldPtr_2d
  end interface

#ifdef CESMCOUPLED
  logical :: cesm_coupled = .true.
  type(ESMF_GeomType_Flag) :: geomtype = ESMF_GEOMTYPE_MESH
#else
  logical :: cesm_coupled = .false.
  type(ESMF_GeomType_Flag) :: geomtype = ESMF_GEOMTYPE_GRID
#endif

  integer            :: rc,dbrc
  integer            :: import_cnt = 0

!===============================================================================
contains
!===============================================================================

  !> This function has a few purposes: 
  !! (1) it imports surface fluxes using data from the mediator; and 
  !! (2) it can apply restoring in SST and SSS.
  !! See \ref section_ocn_import for a summary of the surface fluxes that are
  !! passed from MCT to MOM6, including fluxes that need to be included in the future.

  subroutine mom_import(ocean_public, ocean_grid, importState, ice_ocean_boundary, runtype, rc)

    ! Input/output variables
    type(ocean_public_type)       , intent(in)    :: ocean_public       !< Ocean surface state
    type(ocean_grid_type)         , intent(in)    :: ocean_grid         !< Ocean model grid
    type(ESMF_State)              , intent(inout) :: importState        !< incoming data from mediator
    type(ice_ocean_boundary_type) , intent(inout) :: ice_ocean_boundary !< Ocean boundary forcing
    character(len=*), optional    , intent(in)    :: runtype            !< For cesm only, type of run
    integer                       , intent(inout) :: rc

    ! Local Variables
    type(ESMF_StateItem_Flag)       :: itemFlag
    integer                         :: i, j, n
    integer                         :: isc, iec, jsc, jec
    logical                         :: do_import
    logical                         :: isPresent_lwup
    logical                         :: isPresent_lwdn
    logical                         :: isPresent_lwnet
    character(len=128)              :: fldname
    character(len=128)              :: fldname_x
    character(len=128)              :: fldname_y
    real(ESMF_KIND_R8), allocatable :: taux(:,:)
    real(ESMF_KIND_R8), allocatable :: tauy(:,:)
    character(len=*)  , parameter   :: subname = '(mom_import_cesm)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! -------
    ! import_cnt is used to skip using the import state at the first count for cesm
    ! -------

    if (present(runtype)) then
       import_cnt = import_cnt + 1
       if ((trim(runtype) == 'initial' .and. import_cnt <= 2)) then
          do_import = .false. ! This will skip the first time import information is given
       else
          do_import = .true.
       end if
    else
       do_import = .true.
    end if

    if (do_import) then
       call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

       !----
       ! surface height pressure
       !----
       if (cesm_coupled) then
          fldname = 'Sa_pslv'
       else
          fldname = 'inst_pres_height_surface'
       end if
       call state_getimport(importState, trim(fldname), &
            isc, iec, jsc, jec, ice_ocean_boundary%p, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! near-IR, direct shortwave  (W/m2)
       !----
       call state_getimport(importState, 'mean_net_sw_ir_dir_flx', &
            isc, iec, jsc, jec, ice_ocean_boundary%sw_flux_nir_dir, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! near-IR, diffuse shortwave  (W/m2)
       !----
       call state_getimport(importState, 'mean_net_sw_ir_dif_flx', &
            isc, iec, jsc, jec, ice_ocean_boundary%sw_flux_nir_dif, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! visible, direct shortwave  (W/m2)
       !----
       call state_getimport(importState, 'mean_net_sw_vis_dir_flx', &
            isc, iec, jsc, jec, ice_ocean_boundary%sw_flux_vis_dir, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! visible, diffuse shortwave (W/m2)
       !----
       call state_getimport(importState, 'mean_net_sw_vis_dif_flx', &
            isc, iec, jsc, jec, ice_ocean_boundary%sw_flux_vis_dif, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       ! -------
       ! Net longwave radiation (W/m2)
       ! -------
       ! Different treatment of long wave dependent on atmosphere 
       ! When running with cam or datm - need Foxx_lwup and Faxa_lwdn
       ! When running with fv3 - need mean_net_lw_flx

       call ESMF_StateGet(importState, 'Foxx_lwup', itemFlag, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
       if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
          isPresent_lwup = .true.
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
       else
          isPresent_lwnet = .false.
       end if

       if (isPresent_lwup .and. isPresent_lwdn) then
          ! longwave radiation, sum up and down (W/m2)
          call state_getimport(importState, 'Foxx_lwup',  &
               isc, iec, jsc, jec, ice_ocean_boundary%lw_flux, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out
          call state_getimport(importState, 'Faxa_lwdn', &
               isc, iec, jsc, jec, ice_ocean_boundary%lw_flux, do_sum=.true., rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out
       else if (isPresent_lwnet) then
          ! net longwave radiation, sum up and down (W/m2)
          call state_getimport(importState, 'mean_net_lw_flx',  &
               isc, iec, jsc, jec, ice_ocean_boundary%lw_flux, do_sum=.true., rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out
       end if

       !----
       ! zonal and meridional surface stress
       !----
       if (cesm_coupled) then
          fldname_x = 'Foxx_taux'
          fldname_y = 'Foxx_tauy'
       else
          fldname_x = 'mean_zonal_moment_flx'
          fldname_y = 'mean_merid_moment_flx'
       end if

       allocate (taux(isc:iec,jsc:jec))
       allocate (tauy(isc:iec,jsc:jec))
       call state_getimport(importState, trim(fldname_x), isc, iec, jsc, jec, taux, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
       call state_getimport(importState, trim(fldname_y), isc, iec, jsc, jec, tauy, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       ! rotate taux and tauy from true zonal/meridional to local coordinates
       ! Note - this is the latest calculation from Gustavo - pointed out that the NEMS calculation is incorrect
       if (cesm_coupled) then
          do j = jsc, jec
             do i = isc, iec
                ! TODO (mvertens, 2018-12-28): create a new baseline with these changes
                !ice_ocean_boundary%u_flux(i,j) = ocean_grid%cos_rot(i,j) * taux(i,j) + ocean_grid%sin_rot(i,j) * tauy(i,j)
                !ice_ocean_boundary%v_flux(i,j) = ocean_grid%cos_rot(i,j) * tauy(i,j) - ocean_grid%sin_rot(i,j) * taux(i,j)
                ice_ocean_boundary%u_flux(i,j) = taux(i,j)
                ice_ocean_boundary%v_flux(i,j) = tauy(i,j)
             end do
          end do
       else
          do j = jsc, jec
             do i = isc, iec
                ice_ocean_boundary%u_flux(i,j) = ocean_grid%cos_rot(i,j)*taux(i,j) - ocean_grid%sin_rot(i,j)*tauy(i,j)
                ice_ocean_boundary%v_flux(i,j) = ocean_grid%cos_rot(i,j)*tauy(i,j) + ocean_grid%sin_rot(1,j)*taux(i,j)
             end do
          end do
       end if

       !----
       ! sensible heat flux (W/m2)
       !----
       if (cesm_coupled) then
          fldname = 'Foxx_sen'
       else
          fldname = 'mean_sensi_heat_flx'
       end if
       call state_getimport(importState, trim(fldname), &
            isc, iec, jsc, jec, ice_ocean_boundary%t_flux, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! latent heat flux (W/m2)
       !----
       if (cesm_coupled) then
          ! Note - this field is not exported by the nems mediator
          fldname = 'Foxx_lat'
          call state_getimport(importState, trim(fldname), &
               isc, iec, jsc, jec, ice_ocean_boundary%latent_flux, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out
       end if

       !----
       ! specific humidity flux (W/m2)
       !----
       if (cesm_coupled) then
          fldname = 'Foxx_evap'
       else
          fldname = 'mean_evap_rate'
       end if
       call state_getimport(importState, trim(fldname), &
            isc, iec, jsc, jec, ice_ocean_boundary%q_flux, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! liquid precipitation (rain)
       !----
       if (cesm_coupled) then
          fldname = 'Faxa_rain'
       else
          fldname = 'mean_prec_rate'
       end if
       call state_getimport(importState, trim(fldname), &
            isc, iec, jsc, jec, ice_ocean_boundary%lprec, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! frozen precipitation (snow)
       !----
       if (cesm_coupled) then
          fldname = 'Faxa_snow'
       else
          fldname = 'mean_fprec_rate'
       end if
       call state_getimport(importState, trim(fldname), &
            isc, iec, jsc, jec, ice_ocean_boundary%fprec, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! runoff and heat content of runoff
       !----
       if (cesm_coupled) then
          ! liquid runoff
          fldname = 'Foxx_rofl'
          call state_getimport(importState, trim(fldname),  &
               isc, iec, jsc, jec, ice_ocean_boundary%rofl_flux,rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out

          ! ice runoff
          fldname = 'Foxx_rofi'
          call state_getimport(importState, trim(fldname),  &
               isc, iec, jsc, jec, ice_ocean_boundary%rofi_flux,rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out

          ! GMM, cime does not not have an equivalent for heat_content_lrunoff and
          ! heat_content_frunoff. Setting these to zero for now.
          ice_ocean_boundary%runoff_hflx(:,:)  = 0._ESMF_KIND_R8
          ice_ocean_boundary%calving_hflx(:,:) = 0._ESMF_KIND_R8

       else
          ! total runoff
          fldname = 'mean_runoff_rate'
          call state_getimport(importState, trim(fldname),  &
               isc, iec, jsc, jec, ice_ocean_boundary%runoff, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out

          ! heat content of runoff
          fldname = 'mean_runoff_heat_flux'
          call state_getimport(importState, trim(fldname),  &
               isc, iec, jsc, jec, ice_ocean_boundary%runoff_hflx, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out
       end if

       !----
       ! calving rate and heat flux
       !----
       if (.not. cesm_coupled) then
          fldname = 'mean_calving_rate'
          call state_getimport(importState, trim(fldname),  &
               isc, iec, jsc, jec, ice_ocean_boundary%calving, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out

          fldname = 'mean_calving_heat_flux'
          call state_getimport(importState, trim(fldname),  &
               isc, iec, jsc, jec, ice_ocean_boundary%calving_hflx, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out
       end if

       !----
       ! salt flux from ice
       !----
       if (cesm_coupled) then
          fldname = 'Fioi_salt'
       else
          fldname = 'mean_salt_rate'
       end if
       call state_getimport(importState, trim(fldname),  &
            isc, iec, jsc, jec, ice_ocean_boundary%salt_flux,rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out
       if (cesm_coupled) then
          ! salt flux (minus sign needed here -GMM)  
          ! TODO (mvertens, 2018-12-28): NEMS does not have a minus sign - which one is right?
          do j = jsc,jec
             do i = isc,iec
                ice_ocean_boundary%salt_flux(i,j) = - ice_ocean_boundary%salt_flux(i,j)
             enddo
          enddo
       end if
          
       !----
       ! mass of overlying ice
       !----
       if (.not. cesm_coupled) then
          fldname = 'mass_of_overlying_ice'
          call state_getimport(importState, trim(fldname),  &
               isc, iec, jsc, jec, ice_ocean_boundary%mi, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out
       end if

    end if

  end subroutine mom_import

!===============================================================================

  !> Maps outgoing ocean data to ESMF State
  subroutine mom_export_cesm(ocean_public, ocean_grid, exportState, clock, rc)

    ! Input/output variables
    type(ocean_public_type) , intent(in)    :: ocean_public !< Ocean surface state
    type(ocean_grid_type)   , intent(in)    :: ocean_grid   !< Ocean model grid
    type(ESMF_State)        , intent(inout) :: exportState  !< outgoing data
    type(ESMF_Clock)        , intent(in)    :: clock
    integer                 , intent(inout) :: rc

    ! Local variables
    real            :: ssh(ocean_grid%isd:ocean_grid%ied, ocean_grid%jsd:ocean_grid%jed) !< Local copy of sea_lev with updated halo
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

    !----------------

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
    do j=ocean_grid%jsc, ocean_grid%jec
       jg = j + ocean_grid%jdg_offset
       do i=ocean_grid%isc,ocean_grid%iec
          ig = i + ocean_grid%idg_offset
          n = n+1
          dataPtr_omask(n)   = ocean_grid%mask2dT(i,j)
          dataPtr_t(n)       = ocean_public%t_surf(ig,jg) * ocean_grid%mask2dT(i,j) ! surface temp is in K
          dataPtr_s(n)       = ocean_public%s_surf(ig,jg) * ocean_grid%mask2dT(i,j)
          dataPtr_u(n)       = ocean_public%u_surf(ig,jg) * ocean_grid%mask2dT(i,j)
          dataPtr_v(n)       = ocean_public%v_surf(ig,jg) * ocean_grid%mask2dT(i,j)
          dataPtr_bldepth(n) = ocean_public%OBLD(ig,jg)   * ocean_grid%mask2dT(i,j)
          ! ocean melt and freeze potential (o2x_Fioo_q), W m-2
          if (ocean_public%frazil(ig,jg) > 0.0) then
             ! Frazil: change from J/m^2 to W/m^2
             dataPtr_Fioo_q(n) = ocean_public%frazil(ig,jg) * ocean_grid%mask2dT(i,j) * I_time_int
          else
             ! Melt_potential: change from J/m^2 to W/m^2
             dataPtr_Fioo_q(n) = -ocean_public%melt_potential(ig,jg) * ocean_grid%mask2dT(i,j) * I_time_int !* ncouple_per_day

             ! make sure Melt_potential is always <= 0
             if (dataPtr_Fioo_q(n) > 0.0) dataPtr_Fioo_q(n) = 0.0
          end if
       end do
    end do

    ! Make a copy of ssh in order to do a halo update.
    ! ssh has global indexing with halos

    do j = ocean_grid%jsc, ocean_grid%jec
      jg = j + ocean_grid%jdg_offset
      do i = ocean_grid%isc,ocean_grid%iec
        ig = i + ocean_grid%idg_offset
        ssh(i,j) = ocean_public%sea_lev(ig,jg)
      end do
    end do

    ! Update halo of ssh so we can calculate gradients
    call pass_var(ssh, ocean_grid%domain)

    ! d/dx ssh
    n = 0
    do j=ocean_grid%jsc, ocean_grid%jec
       do i=ocean_grid%isc,ocean_grid%iec
          n = n+1
          ! This is a simple second-order difference
          ! dataPtr_dhdx(n) = 0.5 * (ssh(i+1,j) - ssh(i-1,j)) * ocean_grid%IdxT(i,j) * ocean_grid%mask2dT(i,j)
          ! This is a PLM slope which might be less prone to the A-ocean_grid null mode
          slp_L = (ssh(I,j) - ssh(I-1,j)) * ocean_grid%mask2dCu(I-1,j)
          if (ocean_grid%mask2dCu(I-1,j)==0.) slp_L = 0.
          slp_R = (ssh(I+1,j) - ssh(I,j)) * ocean_grid%mask2dCu(I,j)
          if (ocean_grid%mask2dCu(I+1,j)==0.) slp_R = 0.
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
          dataPtr_dhdx(n) = slope * ocean_grid%IdxT(i,j) * ocean_grid%mask2dT(i,j)
          if (ocean_grid%mask2dT(i,j)==0.) dataPtr_dhdx(n) = 0.0
       enddo
    enddo

    ! d/dy ssh
    n = 0
    do j=ocean_grid%jsc, ocean_grid%jec
       do i=ocean_grid%isc,ocean_grid%iec
          n = n+1
          ! This is a simple second-order difference
          ! dataPtr_dhdy(n) = 0.5 * (ssh(i,j+1) - ssh(i,j-1)) * ocean_grid%IdyT(i,j) * ocean_grid%mask2dT(i,j)
          ! This is a PLM slope which might be less prone to the A-ocean_grid null mode
          slp_L = ssh(i,J) - ssh(i,J-1) * ocean_grid%mask2dCv(i,J-1)
          if (ocean_grid%mask2dCv(i,J-1)==0.) slp_L = 0.

          slp_R = ssh(i,J+1) - ssh(i,J) * ocean_grid%mask2dCv(i,J)
          if (ocean_grid%mask2dCv(i,J+1)==0.) slp_R = 0.

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
          dataPtr_dhdy(n) = slope * ocean_grid%IdyT(i,j) * ocean_grid%mask2dT(i,j)
          if (ocean_grid%mask2dT(i,j)==0.) dataPtr_dhdy(n) = 0.0
       enddo
    enddo

  end subroutine mom_export_cesm

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

    lbnd1 = lbound(dataPtr_dhdx,1)
    ubnd1 = ubound(dataPtr_dhdx,1)
    lbnd2 = lbound(dataPtr_dhdx,2)
    ubnd2 = ubound(dataPtr_dhdx,2)

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

  !===============================================================================

  subroutine State_GetImport(state, fldname, isc, iec, jsc, jec, output, do_sum, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    character(len=*)    , intent(in)    :: fldname
    integer             , intent(in)    :: isc
    integer             , intent(in)    :: iec
    integer             , intent(in)    :: jsc
    integer             , intent(in)    :: jec
    real (ESMF_KIND_R8) , intent(inout) :: output(isc:iec,jsc:jec)
    logical, optional   , intent(in)    :: do_sum
    integer             , intent(out)   :: rc

    ! local variables
    integer                       :: n, i, j, i1, j1
    integer                       :: lbnd1,lbnd2
    real(ESMF_KIND_R8), pointer   :: dataPtr1d(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr2d(:,:)
    character(len=*)  , parameter :: subname='(mom_cap_methods:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (geomtype == ESMF_GEOMTYPE_MESH) then

       ! get field pointer
       call state_getfldptr(state, trim(fldname), dataptr1d, rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       ! determine output array
       n = 0
       do j = jsc,jec
          do i = isc,iec
             n = n + 1 
             if (present(do_sum)) then
                output(i,j)  = output(i,j) + dataPtr1d(n)
             else
                output(i,j)  = dataPtr1d(n)
             end if
          end do
       end do

    else if (geomtype == ESMF_GEOMTYPE_GRID) then

       call state_getfldptr(state, trim(fldname), dataptr2d, rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       lbnd1 = lbound(dataPtr2d,1)
       lbnd2 = lbound(dataPtr2d,2)

       do j = jsc, jec
          j1 = j + lbnd2 - jsc
          do i = isc, iec
             i1 = i + lbnd1 - isc
             if (present(do_sum)) then
                output(i,j) = output(i,j) + dataPtr2d(i1,j1)
             else
                output(i,j) = dataPtr2d(i1,j1)
             end if
          end do
       end do

    end if

  end subroutine State_GetImport

  !===============================================================================

  subroutine State_SetExport(state, fldname, isc, iec, jsc, jec, input, ocean_grid, rc)

    ! ----------------------------------------------
    ! Map input array to export state
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)      , intent(inout) :: state
    character(len=*)      , intent(in)    :: fldname
    integer               , intent(in)    :: isc
    integer               , intent(in)    :: iec
    integer               , intent(in)    :: jsc
    integer               , intent(in)    :: jec
    real (ESMF_KIND_R8)   , intent(in)    :: input(isc:iec,jsc:jec)
    type(ocean_grid_type) , intent(in)    :: ocean_grid
    integer               , intent(out)   :: rc

    ! local variables
    integer                       :: n, i, j, i1, j1, ig,jg
    integer                       :: lbnd1,lbnd2
    real(ESMF_KIND_R8), pointer   :: dataPtr1d(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr2d(:,:)
    character(len=*)  , parameter :: subname='(mom_cap_methods_:state_setimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Indexing notes:
    ! input array from "ocean_public" uses local indexing without halos
    ! mask from "ocean_grid" uses global indexing with halos

    if (geomtype == ESMF_GEOMTYPE_MESH) then

       call state_getfldptr(state, trim(fldname), dataptr1d, rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       n = 0
       do j = jsc, jec
          jg = j + ocean_grid%jsc - jsc
          do i = isc, iec
             ig = i + ocean_grid%isc - isc
             n = n+1
             dataPtr1d(n) = input(i,j) * ocean_grid%mask2dT(ig,jg)
          end do
       end do

    else if (geomtype == ESMF_GEOMTYPE_GRID) then

       call state_getfldptr(state, trim(fldname), dataptr2d, rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       lbnd1 = lbound(dataPtr2d,1)
       lbnd2 = lbound(dataPtr2d,2)

       do j = jsc, jec
          j1 = j + lbnd2 - jsc  
          jg = j + ocean_grid%jsc - jsc
          do i = isc, iec
             i1 = i + lbnd1 - isc
             ig = i + ocean_grid%isc - isc
             dataPtr2d(i1,j1)  = input(i,j) * ocean_grid%mask2dT(ig,jg) 
          end do
       end do

    end if

  end subroutine State_SetExport

end module mom_cap_methods
