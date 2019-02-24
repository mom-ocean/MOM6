module mom_cap_methods

  ! Cap import/export methods for both NEMS and CMEPS

  use ESMF,                only: ESMF_Clock, ESMF_ClockGet, ESMF_time, ESMF_TimeGet
  use ESMF,                only: ESMF_TimeInterval, ESMF_TimeIntervalGet
  use ESMF,                only: ESMF_State, ESMF_StateGet
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
  !DLW
  use MOM_ocean_model,     only: ocean_public_type, ocean_state_type
  use MOM_surface_forcing, only: ice_ocean_boundary_type
  use MOM_grid,            only: ocean_grid_type
  use MOM_domains,         only: pass_var
  use mpp_domains_mod,     only: mpp_get_compute_domain

  ! By default make data private
  implicit none
  private

  ! Public member functions
  public :: mom_set_geomtype
  public :: mom_import
  public :: mom_export

  private :: State_getImport
  private :: State_setExport

  interface State_GetFldPtr
     module procedure State_GetFldPtr_1d
     module procedure State_GetFldPtr_2d
  end interface

  integer                  :: import_cnt = 0
  type(ESMF_GeomType_Flag) :: geomtype

!===============================================================================
contains
!===============================================================================

  subroutine mom_set_geomtype(geomtype_in)
    ! Set module variable geomtype

    type(ESMF_GeomType_Flag), intent(in)    :: geomtype_in     !< mesh or grid

    geomtype = geomtype_in

  end subroutine mom_set_geomtype

!===============================================================================

  !> This function has a few purposes:
  !! (1) it imports surface fluxes using data from the mediator; and
  !! (2) it can apply restoring in SST and SSS.

  subroutine mom_import(ocean_public, ocean_grid, importState, ice_ocean_boundary, runtype, rc)

    ! Input/output variables
    type(ocean_public_type)       , intent(in)    :: ocean_public       !< Ocean surface state
    type(ocean_grid_type)         , intent(in)    :: ocean_grid         !< Ocean model grid
    type(ESMF_State)              , intent(inout) :: importState        !< incoming data from mediator
    type(ice_ocean_boundary_type) , intent(inout) :: ice_ocean_boundary !< Ocean boundary forcing
    character(len=*), optional    , intent(in)    :: runtype            !< For cesm only, type of run
    integer                       , intent(inout) :: rc

    ! Local Variables
    integer                         :: i, j, ig, jg, n
    integer                         :: isc, iec, jsc, jec
    logical                         :: do_import
    character(len=128)              :: fldname
    real(ESMF_KIND_R8), allocatable :: taux(:,:)
    real(ESMF_KIND_R8), allocatable :: tauy(:,:)
    character(len=*)  , parameter   :: subname = '(mom_import)'

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
       ! The following are global indices without halos
       call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

       !----
       ! surface height pressure
       !----
       call state_getimport(importState, 'inst_pres_height_surface', &
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
       call state_getimport(importState, 'mean_net_lw_flx',  &
            isc, iec, jsc, jec, ice_ocean_boundary%lw_flux, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! zonal and meridional surface stress
       !----
       allocate (taux(isc:iec,jsc:jec))
       allocate (tauy(isc:iec,jsc:jec))

       call state_getimport(importState, 'mean_zonal_moment_flx', isc, iec, jsc, jec, taux, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
       call state_getimport(importState, 'mean_merid_moment_flx', isc, iec, jsc, jec, tauy, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       ! rotate taux and tauy from true zonal/meridional to local coordinates
       do j = jsc, jec
          jg = j + ocean_grid%jsc - jsc
          do i = isc, iec
             ig = i + ocean_grid%isc - isc
             ice_ocean_boundary%u_flux(i,j) = ocean_grid%cos_rot(ig,jg)*taux(i,j) &
                                            - ocean_grid%sin_rot(ig,jg)*tauy(i,j)
             ice_ocean_boundary%v_flux(i,j) = ocean_grid%cos_rot(ig,jg)*tauy(i,j) &
                                            + ocean_grid%sin_rot(ig,jg)*taux(i,j)
          end do
       end do

       deallocate(taux, tauy)

       !----
       ! sensible heat flux (W/m2)
       !----
       call state_getimport(importState, 'mean_sensi_heat_flx', &
            isc, iec, jsc, jec, ice_ocean_boundary%t_flux, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! evaporation flux (W/m2)
       !----
       call state_getimport(importState, 'mean_evap_rate', &
            isc, iec, jsc, jec, ice_ocean_boundary%q_flux, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! liquid precipitation (rain)
       !----
       call state_getimport(importState, 'mean_prec_rate', &
            isc, iec, jsc, jec, ice_ocean_boundary%lprec, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! frozen precipitation (snow)
       !----
       call state_getimport(importState, 'mean_fprec_rate', &
            isc, iec, jsc, jec, ice_ocean_boundary%fprec, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! runoff and heat content of runoff
       !----
       ! Note - preset values to 0, if field does not exist in importState, then will simply return
       ! and preset value will be used

       ! liquid runoff
       ice_ocean_boundary%rofl_flux (:,:) = 0._ESMF_KIND_R8
       call state_getimport(importState, 'Foxx_rofl',  &
            isc, iec, jsc, jec, ice_ocean_boundary%rofl_flux,rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       ! ice runoff
       ice_ocean_boundary%rofi_flux (:,:) = 0._ESMF_KIND_R8
       call state_getimport(importState, 'Foxx_rofi',  &
            isc, iec, jsc, jec, ice_ocean_boundary%rofi_flux,rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       ! total runoff
       ice_ocean_boundary%runoff (:,:) = 0._ESMF_KIND_R8
       call state_getimport(importState, 'mean_runoff_rate',  &
            isc, iec, jsc, jec, ice_ocean_boundary%runoff, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       ! heat content of runoff
       ice_ocean_boundary%runoff_hflx(:,:) = 0._ESMF_KIND_R8
       call state_getimport(importState, 'mean_runoff_heat_flux',  &
            isc, iec, jsc, jec, ice_ocean_boundary%runoff_hflx, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! calving rate and heat flux
       !----
       ! Note - preset values to 0, if field does not exist in importState, then will simply return
       ! and preset value will be used

       ice_ocean_boundary%calving(:,:) = 0._ESMF_KIND_R8
       call state_getimport(importState, 'mean_calving_rate',  &
            isc, iec, jsc, jec, ice_ocean_boundary%calving, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       ice_ocean_boundary%calving_hflx(:,:) = 0._ESMF_KIND_R8
       call state_getimport(importState, 'mean_calving_heat_flux',  &
            isc, iec, jsc, jec, ice_ocean_boundary%calving_hflx, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       !----
       ! salt flux from ice
       !----
       call state_getimport(importState, 'mean_salt_rate',  &
            isc, iec, jsc, jec, ice_ocean_boundary%salt_flux,rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       ! TODO: salt flux (minus sign needed here -GMM) - this does not match either NEMS or MCT - so not put in below
       do j = jsc,jec
          do i = isc,iec
             ice_ocean_boundary%salt_flux(i,j) = ice_ocean_boundary%salt_flux(i,j)
          enddo
       enddo

       !----
       ! mass of overlying ice
       !----
       ! Note - preset values to 0, if field does not exist in importState, then will simply return
       ! and preset value will be used

       ice_ocean_boundary%mi(:,:) = 0._ESMF_KIND_R8
       call state_getimport(importState, 'mass_of_overlying_ice',  &
            isc, iec, jsc, jec, ice_ocean_boundary%mi, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

    end if

  end subroutine mom_import

!===============================================================================

  !> Maps outgoing ocean data to ESMF State
  subroutine mom_export(ocean_public, ocean_grid, ocean_state, exportState, clock, rc)

    ! Input/output variables
    type(ocean_public_type) , intent(in)    :: ocean_public !< Ocean surface state
    type(ocean_grid_type)   , intent(in)    :: ocean_grid   !< Ocean model grid
    type(ocean_state_type)  , pointer       :: ocean_state
    type(ESMF_State)        , intent(inout) :: exportState  !< outgoing data
    type(ESMF_Clock)        , intent(in)    :: clock
    integer                 , intent(inout) :: rc

    ! Local variables
    integer                         :: i, j, ig, jg         ! indices
    integer                         :: isc, iec, jsc, jec   ! indices
    integer                         :: iloc, jloc           ! indices
    integer                         :: iglob, jglob         ! indices
    integer                         :: n
    integer                         :: icount
    real                            :: slp_L, slp_R, slp_C
    real                            :: slope, u_min, u_max
    integer                         :: day, secs
    type(ESMF_TimeInterval)         :: timeStep
    integer                         :: dt_int
    real                            :: inv_dt_int  !< The inverse of coupling time interval in s-1.
    type(ESMF_StateItem_Flag)       :: itemFlag
    real(ESMF_KIND_R8), allocatable :: omask(:,:)
    real(ESMF_KIND_R8), allocatable :: melt_potential(:,:)
    real(ESMF_KIND_R8), allocatable :: ocz(:,:), ocm(:,:)
    real(ESMF_KIND_R8), allocatable :: ocz_rot(:,:), ocm_rot(:,:)
    real(ESMF_KIND_R8), allocatable :: ssh(:,:)
    real(ESMF_KIND_R8), allocatable :: dhdx(:,:), dhdy(:,:)
    real(ESMF_KIND_R8), allocatable :: dhdx_rot(:,:), dhdy_rot(:,:)
    character(len=*)  , parameter   :: subname = '(mom_export)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

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
       inv_dt_int = 1.0 / real(dt_int)
    else
       inv_dt_int = 0.0
    end if

    !----------------
    ! Copy from ocean_public to exportstate.
    !----------------

    call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

    ! -------
    ! ocean mask
    ! -------

    !DLW: Retrieve omask from ocean_grid, as for other grid variables
    allocate(omask(isc:iec, jsc:jec))
    do j = jsc, jec
       jg = j + ocean_grid%jsc - jsc
       do i = isc, iec
          ig = i + ocean_grid%isc - isc
          omask(i,j) = nint(ocean_grid%mask2dT(ig,jg))
       enddo
    enddo

    call State_SetExport(exportState, 'ocean_mask', &
         isc, iec, jsc, jec, omask, ocean_grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    deallocate(omask)

    ! -------
    ! Sea surface temperature
    ! -------
    call State_SetExport(exportState, 'sea_surface_temperature', &
         isc, iec, jsc, jec, ocean_public%t_surf, ocean_grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    ! -------
    ! Sea surface salinity
    ! -------
    call State_SetExport(exportState, 's_surf', &
         isc, iec, jsc, jec, ocean_public%s_surf, ocean_grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    ! -------
    ! zonal and meridional currents
    ! -------

    ! rotate ocn current from tripolar grid back to lat/lon grid x,y => latlon (CCW)
    ! "ocean_grid%isc" has no halos and uses local indexing.

    allocate(ocz(isc:iec, jsc:jec))
    allocate(ocm(isc:iec, jsc:jec))
    allocate(ocz_rot(isc:iec, jsc:jec))
    allocate(ocm_rot(isc:iec, jsc:jec))

    do j = jsc, jec
       jg = j + ocean_grid%jsc - jsc
       do i = isc, iec
          ig = i + ocean_grid%isc - isc
          ocz(i,j) = ocean_public%u_surf(i,j)
          ocm(i,j) = ocean_public%v_surf(i,j)
          ocz_rot(i,j) = ocean_grid%cos_rot(ig,jg)*ocz(i,j) + ocean_grid%sin_rot(ig,jg)*ocm(i,j)
          ocm_rot(i,j) = ocean_grid%cos_rot(ig,jg)*ocm(i,j) - ocean_grid%sin_rot(ig,jg)*ocz(i,j)
       end do
    end do

    call State_SetExport(exportState, 'ocn_current_zonal', &
         isc, iec, jsc, jec, ocz_rot, ocean_grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    call State_SetExport(exportState, 'ocn_current_merid', &
         isc, iec, jsc, jec, ocm_rot, ocean_grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    deallocate(ocz, ocm, ocz_rot, ocm_rot)

    ! -------
    ! Boundary layer depth
    ! -------
    call ESMF_StateGet(exportState, 'So_bldepth', itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
       call State_SetExport(exportState, 'So_bldepth', &
            isc, iec, jsc, jec, ocean_public%obld, ocean_grid, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
    end if

    ! -------
    ! Freezing melting potential
    ! -------
    ! melt_potential, defined positive for T>Tfreeze, so need to change sign
    ! Convert from J/m^2 to W/m^2 and make sure Melt_potential is always <= 0

    allocate(melt_potential(isc:iec, jsc:jec))

    do j = jsc,jec
       do i = isc,iec
          if (ocean_public%frazil(i,j) > 0.0) then
             melt_potential(i,j) =  ocean_public%frazil(i,j) * inv_dt_int
          else
             melt_potential(i,j) = -ocean_public%melt_potential(i,j) * inv_dt_int
             if (melt_potential(i,j) > 0.0) melt_potential(i,j) = 0.0
          end if
       end do
    end do

    call State_SetExport(exportState, 'freezing_melting_potential', &
         isc, iec, jsc, jec, melt_potential, ocean_grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    deallocate(melt_potential)

    ! -------
    ! Sea level
    ! -------
    call ESMF_StateGet(exportState, 'sea_level', itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
       call State_SetExport(exportState, 'sea_level', &
            isc, iec, jsc, jec, ocean_public%sea_lev, ocean_grid, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
    end if

    !----------------
    ! Sea-surface zonal and meridional slopes
    !----------------

    allocate(ssh(ocean_grid%isd:ocean_grid%ied,ocean_grid%jsd:ocean_grid%jed)) ! local indices with halos
    allocate(dhdx(isc:iec, jsc:jec))     !global indices without halos
    allocate(dhdy(isc:iec, jsc:jec))     !global indices without halos
    allocate(dhdx_rot(isc:iec, jsc:jec)) !global indices without halos
    allocate(dhdy_rot(isc:iec, jsc:jec)) !global indices without halos

    ssh  = 0.0_ESMF_KIND_R8
    dhdx = 0.0_ESMF_KIND_R8
    dhdy = 0.0_ESMF_KIND_R8

    ! Make a copy of ssh in order to do a halo update (ssh has local indexing with halos)
    do j = ocean_grid%jsc, ocean_grid%jec
       jloc = j + ocean_grid%jdg_offset
       do i = ocean_grid%isc,ocean_grid%iec
          iloc = i + ocean_grid%idg_offset
          ssh(i,j) = ocean_public%sea_lev(iloc,jloc)
       end do
    end do

    ! Update halo of ssh so we can calculate gradients (local indexing)
    call pass_var(ssh, ocean_grid%domain)

    ! d/dx ssh
    ! This is a simple second-order difference
    ! dhdx(i,j) = 0.5 * (ssh(i+1,j) - ssh(i-1,j)) * ocean_grid%IdxT(i,j) * ocean_grid%mask2dT(ig,jg)

    do jglob = jsc, jec
      j  = jglob + ocean_grid%jsc - jsc
      do iglob = isc,iec
        i  = iglob + ocean_grid%isc - isc
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
        dhdx(iglob,jglob) = slope * ocean_grid%IdxT(i,j) * ocean_grid%mask2dT(i,j)
        if (ocean_grid%mask2dT(i,j)==0.) dhdx(iglob,jglob) = 0.0
      end do
    end do

    ! d/dy ssh
    ! This is a simple second-order difference
    ! dhdy(i,j) = 0.5 * (ssh(i,j+1) - ssh(i,j-1)) * ocean_grid%IdyT(i,j) * ocean_grid%mask2dT(ig,jg)

    do jglob = jsc, jec
      j = jglob + ocean_grid%jsc - jsc
      do iglob = isc,iec
         i = iglob + ocean_grid%isc - isc
        ! This is a PLM slope which might be less prone to the A-ocean_grid null mode
        slp_L = ssh(i,J) - ssh(i,J-1) * ocean_grid%mask2dCv(i,j-1)
        if (ocean_grid%mask2dCv(i,j-1)==0.) slp_L = 0.
        slp_R = ssh(i,J+1) - ssh(i,J) * ocean_grid%mask2dCv(i,j)
        if (ocean_grid%mask2dCv(i,j+1)==0.) slp_R = 0.
        slp_C = 0.5 * (slp_L + slp_R)
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
        dhdy(iglob,jglob) = slope * ocean_grid%IdyT(i,j) * ocean_grid%mask2dT(i,j)
        if (ocean_grid%mask2dT(i,j)==0.) dhdy(iglob,jglob) = 0.0
      end do
    end do

    ! rotate slopes from tripolar grid back to lat/lon grid,  x,y => latlon (CCW)
    ! "ocean_grid" uses has halos and uses global indexing.

    do j = jsc, jec
       jg = j + ocean_grid%jsc - jsc
       do i = isc, iec
          ig = i + ocean_grid%isc - isc
          dhdx_rot(i,j) = ocean_grid%cos_rot(ig,jg)*dhdx(i,j) + ocean_grid%sin_rot(ig,jg)*dhdy(i,j)
          dhdy_rot(i,j) = ocean_grid%cos_rot(ig,jg)*dhdy(i,j) - ocean_grid%sin_rot(ig,jg)*dhdx(i,j)
       end do
    end do

    call State_SetExport(exportState, 'sea_surface_slope_zonal', &
         isc, iec, jsc, jec, dhdx_rot, ocean_grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    call State_SetExport(exportState, 'sea_surface_slope_merid', &
         isc, iec, jsc, jec, dhdy_rot, ocean_grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    deallocate(ssh, dhdx, dhdy, dhdx_rot, dhdy_rot)

  end subroutine mom_export

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
    type(ESMF_StateItem_Flag)     :: itemFlag
    integer                       :: n, i, j, i1, j1
    integer                       :: lbnd1,lbnd2
    real(ESMF_KIND_R8), pointer   :: dataPtr1d(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr2d(:,:)
    character(len=*)  , parameter :: subname='(mom_cap_methods:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, trim(fldname), itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then

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
    type(ESMF_StateItem_Flag)     :: itemFlag
    integer                       :: n, i, j, i1, j1, ig,jg
    integer                       :: lbnd1,lbnd2
    real(ESMF_KIND_R8), pointer   :: dataPtr1d(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr2d(:,:)
    character(len=*)  , parameter :: subname='(mom_cap_methods:state_setexport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Indexing notes:
    ! input array from "ocean_public" uses local indexing without halos
    ! mask from "ocean_grid" uses global indexing with halos

    call ESMF_StateGet(State, trim(fldname), itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then

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

    end if

  end subroutine State_SetExport

end module mom_cap_methods
