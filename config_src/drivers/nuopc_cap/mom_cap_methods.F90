!> Contains import/export methods for CMEPS.
module MOM_cap_methods

use ESMF,                      only: ESMF_Clock, ESMF_ClockGet, ESMF_time, ESMF_TimeGet
use ESMF,                      only: ESMF_TimeInterval, ESMF_TimeIntervalGet
use ESMF,                      only: ESMF_State, ESMF_StateGet
use ESMF,                      only: ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
use ESMF,                      only: ESMF_GridComp, ESMF_Mesh, ESMF_MeshGet, ESMF_Grid, ESMF_GridCreate
use ESMF,                      only: ESMF_DistGrid, ESMF_DistGridCreate
use ESMF,                      only: ESMF_KIND_R8, ESMF_SUCCESS, ESMF_LogFoundError
use ESMF,                      only: ESMF_LOGERR_PASSTHRU, ESMF_LOGMSG_INFO, ESMF_LOGWRITE
use ESMF,                      only: ESMF_LogSetError, ESMF_RC_MEM_ALLOCATE
use ESMF,                      only: ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
use ESMF,                      only: ESMF_GEOMTYPE_FLAG, ESMF_GEOMTYPE_GRID, ESMF_GEOMTYPE_MESH
use ESMF,                      only: ESMF_RC_VAL_OUTOFRANGE, ESMF_INDEX_DELOCAL, ESMF_MESHLOC_ELEMENT
use ESMF,                      only: ESMF_TYPEKIND_R8, ESMF_FIELDSTATUS_COMPLETE
use ESMF,                      only: ESMF_FieldStatus_Flag, ESMF_LOGMSG_ERROR, ESMF_FAILURE, ESMF_MAXSTR
use ESMF,                      only: operator(/=), operator(==)
use MOM_ocean_model_nuopc,     only: ocean_public_type, ocean_state_type
use MOM_surface_forcing_nuopc, only: ice_ocean_boundary_type
use MOM_grid,                  only: ocean_grid_type
use MOM_domains,               only: pass_var
use mpp_domains_mod,           only: mpp_get_compute_domain

! By default make data private
implicit none; private

! Public member functions
public :: mom_set_geomtype
public :: mom_import
public :: mom_export
public :: state_diagnose
public :: ChkErr

interface State_getImport
   module procedure State_getImport_2d
   module procedure State_getImport_3d ! third dimension being an ungridded dimension
end interface

private :: State_setExport

!> Get field pointer
interface State_GetFldPtr
  module procedure State_GetFldPtr_1d
  module procedure State_GetFldPtr_2d
end interface

integer                  :: import_cnt = 0!< used to skip using the import state
                                          !! at the first count for cesm
type(ESMF_GeomType_Flag) :: geomtype      !< SMF type describing type of
                                          !! geometry (mesh or grid)

! area correction factors for fluxes send and received from mediator
! these actors are ONLY valid for meshes that are read in - so do not need them for
! grids that are calculated internally

real(ESMF_KIND_R8), public, allocatable :: mod2med_areacor(:) ! ratios of model areas to input mesh areas
real(ESMF_KIND_R8), public, allocatable :: med2mod_areacor(:) ! ratios of input mesh areas to model areas
character(len=*),parameter :: u_FILE_u =  __FILE__

contains

!> Sets module variable geometry type
subroutine mom_set_geomtype(geomtype_in)
  type(ESMF_GeomType_Flag), intent(in)    :: geomtype_in !< ESMF type describing type of
                                                         !! geometry (mesh or grid)

  geomtype = geomtype_in

end subroutine mom_set_geomtype

!> This function has a few purposes:
!! (1) it imports surface fluxes using data from the mediator; and
!! (2) it can apply restoring in SST and SSS.
subroutine mom_import(ocean_public, ocean_grid, importState, ice_ocean_boundary, rc)
  type(ocean_public_type)       , intent(in)    :: ocean_public       !< Ocean surface state
  type(ocean_grid_type)         , intent(in)    :: ocean_grid         !< Ocean model grid
  type(ESMF_State)              , intent(inout) :: importState        !< incoming data from mediator
  type(ice_ocean_boundary_type) , intent(inout) :: ice_ocean_boundary !< Ocean boundary forcing
  integer                       , intent(inout) :: rc                 !< Return code

  ! Local Variables
  integer                         :: i, j, ib, ig, jg, n
  integer                         :: isc, iec, jsc, jec
  integer                         :: esmf_ind
  integer                         :: nsc ! number of stokes drift components
  character(len=128)              :: fldname
  real(ESMF_KIND_R8), allocatable :: taux(:,:)
  real(ESMF_KIND_R8), allocatable :: tauy(:,:)
  real(ESMF_KIND_R8), allocatable :: stkx(:,:,:)
  real(ESMF_KIND_R8), allocatable :: stky(:,:,:)
  logical                         :: med_has_co2
  character(len=*)  , parameter   :: subname = '(mom_import)'

  rc = ESMF_SUCCESS

  ! -------
  ! import_cnt is used to skip using the import state at the first count for cesm
  ! -------

  ! The following are global indices without halos
  call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

  !----
  ! surface height pressure
  !----
  call state_getimport(importState, 'Sa_pslv', isc, iec, jsc, jec, ice_ocean_boundary%p, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! near-IR, direct shortwave  (W/m2)
  !----
  call state_getimport(importState, 'Foxx_swnet_idr', isc, iec, jsc, jec, &
       ice_ocean_boundary%sw_flux_nir_dir, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! near-IR, diffuse shortwave  (W/m2)
  !----
  call state_getimport(importState, 'Foxx_swnet_idf', isc, iec, jsc, jec, &
       ice_ocean_boundary%sw_flux_nir_dif, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! visible, direct shortwave  (W/m2)
  !----
  call state_getimport(importState, 'Foxx_swnet_vdr', isc, iec, jsc, jec, &
       ice_ocean_boundary%sw_flux_vis_dir, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! visible, diffuse shortwave (W/m2)
  !----
  call state_getimport(importState, 'Foxx_swnet_vdf', isc, iec, jsc, jec, &
       ice_ocean_boundary%sw_flux_vis_dif, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! -------
  ! Net longwave radiation (W/m2)
  ! -------
  call state_getimport(importState, 'Foxx_lwnet', isc, iec, jsc, jec, &
       ice_ocean_boundary%lw_flux, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! zonal and meridional surface stress
  !----
  allocate (taux(isc:iec,jsc:jec))
  allocate (tauy(isc:iec,jsc:jec))

  call state_getimport(importState, 'Foxx_taux', isc, iec, jsc, jec, taux, &
       areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  call state_getimport(importState, 'Foxx_tauy', isc, iec, jsc, jec, tauy, &
       areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! rotate taux and tauy from true zonal/meridional to local coordinates
  do j = jsc, jec
    jg = j + ocean_grid%jsc - jsc
    do i = isc, iec
      ig = i + ocean_grid%isc - isc
      ice_ocean_boundary%u_flux(i,j) = ocean_grid%cos_rot(ig,jg)*taux(i,j) &
           - ocean_grid%sin_rot(ig,jg)*tauy(i,j)
      ice_ocean_boundary%v_flux(i,j) = ocean_grid%cos_rot(ig,jg)*tauy(i,j) &
           + ocean_grid%sin_rot(ig,jg)*taux(i,j)
    enddo
  enddo

  deallocate(taux, tauy)

  !----
  ! sensible heat flux (W/m2)
  !----
  call state_getimport(importState, 'Foxx_sen', isc, iec, jsc, jec, &
       ice_ocean_boundary%t_flux, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! evaporation flux (W/m2)
  !----
  call state_getimport(importState, 'Foxx_evap', isc, iec, jsc, jec, &
       ice_ocean_boundary%q_flux, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! liquid precipitation (rain)
  !----
  call state_getimport(importState, 'Faxa_rain', isc, iec, jsc, jec, &
       ice_ocean_boundary%lprec, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! frozen precipitation (snow)
  !----
  call state_getimport(importState, 'Faxa_snow', isc, iec, jsc, jec, &
       ice_ocean_boundary%fprec, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! mass and heat content of liquid and frozen runoff
  !----
  ! Note - preset values to 0, if field does not exist in importState, then will simply return
  ! and preset value will be used

  ! liquid runoff
  ice_ocean_boundary%lrunoff (:,:) = 0._ESMF_KIND_R8
  call state_getimport(importState, 'Foxx_rofl',  &
       isc, iec, jsc, jec, ice_ocean_boundary%lrunoff, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! ice runoff
  ice_ocean_boundary%frunoff (:,:) = 0._ESMF_KIND_R8
  call state_getimport(importState, 'Foxx_rofi',  &
       isc, iec, jsc, jec, ice_ocean_boundary%frunoff, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! liquid glc runoff
  if ( associated(ice_ocean_boundary%lrunoff_glc) ) then
    ice_ocean_boundary%lrunoff_glc (:,:) = 0._ESMF_KIND_R8
    call state_getimport(importState, 'Forr_rofl_glc',  &
         isc, iec, jsc, jec, ice_ocean_boundary%lrunoff_glc, areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  ! frozen glc runoff
  if ( associated(ice_ocean_boundary%frunoff_glc) ) then
    ice_ocean_boundary%frunoff_glc (:,:) = 0._ESMF_KIND_R8
    call state_getimport(importState, 'Forr_rofi_glc',  &
         isc, iec, jsc, jec, ice_ocean_boundary%frunoff_glc, areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  !----
  ! Enthalpy terms
  !----

  !----
  ! enthalpy from liquid precipitation (hrain)
  !----
  if ( associated(ice_ocean_boundary%hrain) ) then
    call state_getimport(importState, 'Foxx_hrain', isc, iec, jsc, jec, &
         ice_ocean_boundary%hrain, areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  end if

  !----
  ! enthalpy from frozen precipitation (hsnow)
  !----
  if ( associated(ice_ocean_boundary%hsnow) ) then
    call state_getimport(importState, 'Foxx_hsnow', isc, iec, jsc, jec, &
         ice_ocean_boundary%hsnow, areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  end if

  !----
  ! enthalpy from liquid runoff (hrofl)
  !----
  if ( associated(ice_ocean_boundary%hrofl) ) then
    call state_getimport(importState, 'Foxx_hrofl', isc, iec, jsc, jec, &
         ice_ocean_boundary%hrofl, areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  end if

  !----
  ! enthalpy from frozen runoff (hrofi)
  !----
  if ( associated(ice_ocean_boundary%hrofi) ) then
    call state_getimport(importState, 'Foxx_hrofi', isc, iec, jsc, jec, &
         ice_ocean_boundary%hrofi, areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  end if

  !----
  ! enthalpy from liquid glc runoff (hrofl_glc)
  !----
  if ( associated(ice_ocean_boundary%hrofl_glc) ) then
    call state_getimport(importState, 'Foxx_hrofl_glc', isc, iec, jsc, jec, &
         ice_ocean_boundary%hrofl_glc, areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  end if

  !----
  ! enthalpy from frozen glc runoff (hrofi_glc)
  !----
  if ( associated(ice_ocean_boundary%hrofi_glc) ) then
    call state_getimport(importState, 'Foxx_hrofi_glc', isc, iec, jsc, jec, &
         ice_ocean_boundary%hrofi_glc, areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  end if
  !----
  ! enthalpy from evaporation (hevap)
  !----
  if ( associated(ice_ocean_boundary%hevap) ) then
    call state_getimport(importState, 'Foxx_hevap', isc, iec, jsc, jec, &
         ice_ocean_boundary%hevap, areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  end if

  !----
  ! enthalpy from condensation (hcond)
  !----
  if ( associated(ice_ocean_boundary%hcond) ) then
    call state_getimport(importState, 'Foxx_hcond', isc, iec, jsc, jec, &
         ice_ocean_boundary%hcond, areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  !--------------!
  ! MARBL fields !
  !--------------!

  ! seaice_dust_flux, nhx_dep, and noy_dep are single fields from the coupler
  ! atm_fine_dust_flux, atm_coarse_dust_flux, atm_bc_flux, and seaice_bc_flux
  ! are all sums of multiple fields and will be treated slightly differently
  ! For those fields, we use do_sum = .true.

  !----
  ! nhx deposition
  !----
  if (associated(ice_ocean_boundary%nhx_dep)) then
    call state_getimport(importState, 'Faxa_ndep',  &
        isc, iec, jsc, jec, ice_ocean_boundary%nhx_dep(:,:), &
        areacor=med2mod_areacor, esmf_ind=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

    !----
    ! noy deposition
    !----
  if (associated(ice_ocean_boundary%noy_dep)) then
    call state_getimport(importState, 'Faxa_ndep',  &
        isc, iec, jsc, jec, ice_ocean_boundary%noy_dep(:,:), &
        areacor=med2mod_areacor, esmf_ind=2, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  !----
  ! atmospheric CO2 concentration
  ! might not be passed from atmosphere component,
  ! in which the pointer(s) will not be associated
  !----
  if (associated(ice_ocean_boundary%atm_co2_prog)) then
    call state_getimport(importState, 'Sa_co2prog',  &
        isc, iec, jsc, jec, ice_ocean_boundary%atm_co2_prog(:,:), rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif
  if (associated(ice_ocean_boundary%atm_co2_diag)) then
    call state_getimport(importState, 'Sa_co2diag',  &
        isc, iec, jsc, jec, ice_ocean_boundary%atm_co2_diag(:,:), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  endif

  !----
  ! fine dust flux from atmosphere
  !----
  if (associated(ice_ocean_boundary%atm_fine_dust_flux)) then
    ice_ocean_boundary%atm_fine_dust_flux(:,:) = 0._ESMF_KIND_R8
    call state_getimport(importState, 'Faxa_dstwet', &
        isc, iec, jsc, jec, ice_ocean_boundary%atm_fine_dust_flux(:,:), &
        areacor=med2mod_areacor, do_sum=.true., esmf_ind=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Faxa_dstdry',  &
        isc, iec, jsc, jec, ice_ocean_boundary%atm_fine_dust_flux(:,:), &
        areacor=med2mod_areacor, do_sum=.true., esmf_ind=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  !----
  ! coarse dust flux from atmosphere
  !----
  if (associated(ice_ocean_boundary%atm_coarse_dust_flux)) then
    ice_ocean_boundary%atm_coarse_dust_flux(:,:) = 0._ESMF_KIND_R8
    do esmf_ind=2,4
      call state_getimport(importState, 'Faxa_dstwet',  &
          isc, iec, jsc, jec, ice_ocean_boundary%atm_coarse_dust_flux(:,:), &
          areacor=med2mod_areacor, do_sum=.true., esmf_ind=esmf_ind, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call state_getimport(importState, 'Faxa_dstdry',  &
          isc, iec, jsc, jec, ice_ocean_boundary%atm_coarse_dust_flux(:,:), &
          areacor=med2mod_areacor, do_sum=.true., esmf_ind=esmf_ind, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo
  end if

  !----
  ! dust flux from sea ice
  !----
  if (associated(ice_ocean_boundary%seaice_dust_flux)) then
    call state_getimport(importState, 'Fioi_flxdst',  &
        isc, iec, jsc, jec, ice_ocean_boundary%seaice_dust_flux, &
        areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  !----
  ! black carbon flux from atmosphere
  !----
  if (associated(ice_ocean_boundary%atm_bc_flux)) then
    ice_ocean_boundary%atm_bc_flux(:,:) = 0._ESMF_KIND_R8
    do esmf_ind=1,3
      call state_getimport(importState, 'Faxa_bcph',  &
          isc, iec, jsc, jec, ice_ocean_boundary%atm_bc_flux(:,:), &
          areacor=med2mod_areacor, do_sum=.true., esmf_ind=esmf_ind, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo
  endif

  !----
  ! black carbon flux from sea ice
  !----
  if (associated(ice_ocean_boundary%seaice_bc_flux)) then
    ice_ocean_boundary%seaice_bc_flux(:,:) = 0._ESMF_KIND_R8
    call state_getimport(importState, 'Fioi_bcpho',  &
        isc, iec, jsc, jec, ice_ocean_boundary%seaice_bc_flux(:,:), &
        areacor=med2mod_areacor, do_sum=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Fioi_bcphi',  &
        isc, iec, jsc, jec, ice_ocean_boundary%seaice_bc_flux(:,:), &
        areacor=med2mod_areacor, do_sum=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  ! Fields coming from coupler per ice category
  if (ice_ocean_boundary%ice_ncat > 0) then
    call state_getimport(importState, 'Sf_afracr',  &
        isc, iec, jsc, jec, ice_ocean_boundary%afracr(:,:), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    call state_getimport(importState, 'Foxx_swnet_afracr',  &
        isc, iec, jsc, jec, ice_ocean_boundary%swnet_afracr(:,:), &
        areacor=med2mod_areacor, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    call state_getimport(importState, 'Fioi_swpen_ifrac_n',  &
        isc, iec, jsc, jec, 1, ice_ocean_boundary%ice_ncat, &
        ice_ocean_boundary%swpen_ifrac_n(:,:,:), &
        areacor=med2mod_areacor, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    call state_getimport(importState, 'Si_ifrac_n',  &
        isc, iec, jsc, jec, 1, ice_ocean_boundary%ice_ncat, &
        ice_ocean_boundary%ifrac_n(:,:,:), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
  endif ! multiple ice categories

  !----
  ! salt flux from ice
  !----
  ice_ocean_boundary%salt_flux(:,:) = 0._ESMF_KIND_R8
  call state_getimport(importState, 'Fioi_salt', isc, iec, jsc, jec, &
       ice_ocean_boundary%salt_flux, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! snow&ice melt heat flux  (W/m^2)
  !----
  ice_ocean_boundary%seaice_melt_heat(:,:) = 0._ESMF_KIND_R8
  call state_getimport(importState, 'Fioi_melth', isc, iec, jsc, jec, &
       ice_ocean_boundary%seaice_melt_heat, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! snow&ice melt water flux  (W/m^2)
  !----
  ice_ocean_boundary%seaice_melt(:,:) = 0._ESMF_KIND_R8
  call state_getimport(importState, 'Fioi_meltw', isc, iec, jsc, jec, &
       ice_ocean_boundary%seaice_melt, areacor=med2mod_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! mass of overlying ice
  !----
  ! Note - preset values to 0, if field does not exist in importState, then will simply return
  ! and preset value will be used
  ice_ocean_boundary%mi(:,:) = 0._ESMF_KIND_R8
  call state_getimport(importState, 'mass_of_overlying_ice', isc, iec, jsc, jec, &
       ice_ocean_boundary%mi,rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! sea-ice fraction
  !----
  ice_ocean_boundary%ice_fraction(:,:) = 0._ESMF_KIND_R8
  call state_getimport(importState, 'Si_ifrac', isc, iec, jsc, jec, &
       ice_ocean_boundary%ice_fraction, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! 10m wind squared
  !----
  ice_ocean_boundary%u10_sqr(:,:) = 0._ESMF_KIND_R8
  call state_getimport(importState, 'So_duu10n', isc, iec, jsc, jec, &
       ice_ocean_boundary%u10_sqr, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !----
  ! Langmuir enhancement factor
  !----
  if ( associated(ice_ocean_boundary%lamult) ) then
    ice_ocean_boundary%lamult (:,:) = 0._ESMF_KIND_R8
    call state_getimport(importState, 'Sw_lamult', isc, iec, jsc, jec, &
         ice_ocean_boundary%lamult, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  !----
  ! Partitioned Stokes Drift Components
  !----
  if ( associated(ice_ocean_boundary%ustkb) ) then
      nsc = Ice_ocean_boundary%num_stk_bands
      allocate(stkx(isc:iec,jsc:jec,1:nsc))
      allocate(stky(isc:iec,jsc:jec,1:nsc))

      call state_getimport(importState,'Sw_pstokes_x', isc, iec, jsc, jec, 1, nsc, stkx, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call state_getimport(importState,'Sw_pstokes_y', isc, iec, jsc, jec, 1, nsc, stky, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! rotate from true zonal/meridional to local coordinates
      do j = jsc, jec
        jg = j + ocean_grid%jsc - jsc
        do i = isc, iec
          ig = i + ocean_grid%isc - isc
          !rotate
          do ib = 1, nsc
            ice_ocean_boundary%ustkb(i,j,ib) = ocean_grid%cos_rot(ig,jg)*stkx(i,j,ib) &
                 - ocean_grid%sin_rot(ig,jg)*stky(i,j,ib)
            ice_ocean_boundary%vstkb(i,j,ib) = ocean_grid%cos_rot(ig,jg)*stky(i,j,ib) &
                 + ocean_grid%sin_rot(ig,jg)*stkx(i,j,ib)
          enddo
          ! apply masks
          ice_ocean_boundary%ustkb(i,j,:) = ice_ocean_boundary%ustkb(i,j,:) * ocean_grid%mask2dT(ig,jg)
          ice_ocean_boundary%vstkb(i,j,:) = ice_ocean_boundary%vstkb(i,j,:) * ocean_grid%mask2dT(ig,jg)
        enddo
      enddo
      deallocate(stkx,stky)
  endif

end subroutine mom_import

!> Maps outgoing ocean data to ESMF State
subroutine mom_export(ocean_public, ocean_grid, ocean_state, exportState, clock, rc)
  type(ocean_public_type) , intent(in)    :: ocean_public !< Ocean surface state
  type(ocean_grid_type)   , intent(in)    :: ocean_grid   !< Ocean model grid
  type(ocean_state_type)  , pointer       :: ocean_state  !< Ocean state
  type(ESMF_State)        , intent(inout) :: exportState  !< outgoing data
  type(ESMF_Clock)        , intent(in)    :: clock        !< ESMF clock
  integer                 , intent(inout) :: rc           !< Return code

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

  rc = ESMF_SUCCESS

  call ESMF_ClockGet( clock, timeStep=timeStep, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_TimeIntervalGet( timeStep, s=dt_int, rc=rc )
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! Use Adcroft's rule of reciprocals; it does the right thing here.
  if (real(dt_int) > 0.0) then
    inv_dt_int = 1.0 / real(dt_int)
  else
    inv_dt_int = 0.0
  endif

  !----------------
  ! Copy from ocean_public to exportstate.
  !----------------

  call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

  ! -------
  ! ocean mask
  ! -------

  allocate(omask(isc:iec, jsc:jec))
  do j = jsc, jec
    jg = j + ocean_grid%jsc - jsc
    do i = isc, iec
      ig = i + ocean_grid%isc - isc
      omask(i,j) = nint(ocean_grid%mask2dT(ig,jg))
    enddo
  enddo

  call State_SetExport(exportState, 'So_omask', isc, iec, jsc, jec, omask, ocean_grid, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  deallocate(omask)

  ! -------
  ! Sea surface temperature
  ! -------
  call State_SetExport(exportState, 'So_t', isc, iec, jsc, jec, ocean_public%t_surf, ocean_grid, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! -------
  ! Sea surface salinity
  ! -------
  call State_SetExport(exportState, 'So_s', isc, iec, jsc, jec, ocean_public%s_surf, ocean_grid, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! -------
  ! zonal and meridional currents
  ! -------

  ! rotate ocn current from tripolar grid back to lat/lon grid x,y => latlon (CCW)
  ! "ocean_grid" has halos and uses local indexing.

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
    enddo
  enddo

  call State_SetExport(exportState, 'So_u', isc, iec, jsc, jec, ocz_rot, ocean_grid, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call State_SetExport(exportState, 'So_v', isc, iec, jsc, jec, ocm_rot, ocean_grid, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  deallocate(ocz, ocm, ocz_rot, ocm_rot)

  ! -------
  ! Boundary layer depth
  ! -------
  call ESMF_StateGet(exportState, 'So_bldepth', itemFlag, rc=rc)
  if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
    call State_SetExport(exportState, 'So_bldepth', isc, iec, jsc, jec, &
         ocean_public%obld, ocean_grid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

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
      endif
    enddo
  enddo

  call State_SetExport(exportState, 'Fioo_q', isc, iec, jsc, jec, &
       melt_potential, ocean_grid, areacor=mod2med_areacor, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  deallocate(melt_potential)

  ! -------
  ! Sea level
  ! -------
  call ESMF_StateGet(exportState, 'sea_level', itemFlag, rc=rc)
  if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
    call State_SetExport(exportState, 'sea_level', &
         isc, iec, jsc, jec, ocean_public%sea_lev, ocean_grid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  !----------------
  ! Sea-surface zonal and meridional slopes
  !----------------

  allocate(ssh(ocean_grid%isd:ocean_grid%ied,ocean_grid%jsd:ocean_grid%jed), & ! local indices with halos
           dhdx(isc:iec, jsc:jec),                                           & !global indices without halos
           dhdy(isc:iec, jsc:jec),                                           & !global indices without halos
           source=0.0_ESMF_KIND_R8)
  allocate(dhdx_rot(isc:iec, jsc:jec)) !global indices without halos
  allocate(dhdy_rot(isc:iec, jsc:jec)) !global indices without halos

  ! Make a copy of ssh in order to do a halo update (ssh has local indexing with halos)
  do j = ocean_grid%jsc, ocean_grid%jec
    jloc = j + ocean_grid%jdg_offset
    do i = ocean_grid%isc,ocean_grid%iec
      iloc = i + ocean_grid%idg_offset
      ssh(i,j) = ocean_public%sea_lev(iloc,jloc)
    enddo
  enddo

  ! Update halo of ssh so we can calculate gradients (local indexing)
  call pass_var(ssh, ocean_grid%domain)

  ! d/dx ssh
  ! This is a simple second-order difference
  ! dhdx(i,j) = 0.5 * (ssh(i+1,j) - ssh(i-1,j)) * ocean_grid%US%m_to_L*ocean_grid%IdxT(i,j) * ocean_grid%mask2dT(ig,jg)

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
      endif
      dhdx(iglob,jglob) = slope * ocean_grid%US%m_to_L*ocean_grid%IdxT(i,j) * ocean_grid%mask2dT(i,j)
      if (ocean_grid%mask2dT(i,j)==0.) dhdx(iglob,jglob) = 0.0
    enddo
  enddo

  ! d/dy ssh
  ! This is a simple second-order difference
  ! dhdy(i,j) = 0.5 * (ssh(i,j+1) - ssh(i,j-1)) * ocean_grid%US%m_to_L*ocean_grid%IdyT(i,j) * ocean_grid%mask2dT(ig,jg)

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
      endif
      dhdy(iglob,jglob) = slope * ocean_grid%US%m_to_L*ocean_grid%IdyT(i,j) * ocean_grid%mask2dT(i,j)
      if (ocean_grid%mask2dT(i,j)==0.) dhdy(iglob,jglob) = 0.0
    enddo
  enddo

  ! rotate slopes from tripolar grid back to lat/lon grid,  x,y => latlon (CCW)
  ! "ocean_grid" uses has halos and uses local indexing.

  do j = jsc, jec
    jg = j + ocean_grid%jsc - jsc
    do i = isc, iec
      ig = i + ocean_grid%isc - isc
      dhdx_rot(i,j) = ocean_grid%cos_rot(ig,jg)*dhdx(i,j) + ocean_grid%sin_rot(ig,jg)*dhdy(i,j)
      dhdy_rot(i,j) = ocean_grid%cos_rot(ig,jg)*dhdy(i,j) - ocean_grid%sin_rot(ig,jg)*dhdx(i,j)
    enddo
  enddo

  call State_SetExport(exportState, 'So_dhdx', isc, iec, jsc, jec, dhdx_rot, ocean_grid, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call State_SetExport(exportState, 'So_dhdy', isc, iec, jsc, jec, dhdy_rot, ocean_grid, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! -------
  ! CO2 Flux
  ! -------
  call ESMF_StateGet(exportState, 'Faoo_fco2_ocn', itemFlag, rc=rc)
  if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
    call State_SetExport(exportState, 'Faoo_fco2_ocn', isc, iec, jsc, jec, &
        ocean_public%fco2_ocn, ocean_grid, areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  deallocate(ssh, dhdx, dhdy, dhdx_rot, dhdy_rot)

end subroutine mom_export

!> Get field pointer 1D
subroutine State_GetFldPtr_1d(State, fldname, fldptr, rc)
  type(ESMF_State)            , intent(in)  :: State    !< ESMF state
  character(len=*)            , intent(in)  :: fldname  !< Field name
  real(ESMF_KIND_R8), pointer , intent(in)  :: fldptr(:)!< Pointer to the 1D field
  integer, optional           , intent(out) :: rc       !< Return code

  ! local variables
  type(ESMF_Field) :: lfield
  integer :: lrc
  character(len=*),parameter :: subname='(MOM_cap:State_GetFldPtr)'

  call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=lrc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  if (present(rc)) rc = lrc

end subroutine State_GetFldPtr_1d

!> Get field pointer 2D
subroutine State_GetFldPtr_2d(State, fldname, fldptr, rc)
  type(ESMF_State)            , intent(in)  :: State      !< ESMF state
  character(len=*)            , intent(in)  :: fldname    !< Field name
  real(ESMF_KIND_R8), pointer , intent(in)  :: fldptr(:,:)!< Pointer to the 2D field
  integer, optional           , intent(out) :: rc         !< Return code

  ! local variables
  type(ESMF_Field) :: lfield
  integer :: lrc
  character(len=*),parameter :: subname='(MOM_cap:State_GetFldPtr)'

  call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=lrc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  if (present(rc)) rc = lrc

end subroutine State_GetFldPtr_2d

!> Map 2d import state field to output array
subroutine State_GetImport_2d(state, fldname, isc, iec, jsc, jec, output, do_sum, areacor, esmf_ind, rc)
  type(ESMF_State)    , intent(in)    :: state   !< ESMF state
  character(len=*)    , intent(in)    :: fldname !< Field name
  integer             , intent(in)    :: isc     !< The start i-index of cell centers within
                                                 !! the computational domain
  integer             , intent(in)    :: iec     !< The end i-index of cell centers within the
                                                 !! computational domain
  integer             , intent(in)    :: jsc     !< The start j-index of cell centers within
                                                 !! the computational domain
  integer             , intent(in)    :: jec     !< The end j-index of cell centers within
                                                 !! the computational domain
  real (ESMF_KIND_R8) , intent(inout) :: output(isc:iec,jsc:jec)!< Output 2D array
  logical, optional   , intent(in)    :: do_sum  !< If true, sums the data
  real (ESMF_KIND_R8), optional,  intent(in) :: areacor(:) !< flux area correction factors
                                                           !! applicable to meshes
  integer,             optional, intent(in) :: esmf_ind
  integer             , intent(out)   :: rc      !< Return code

  ! local variables
  type(ESMF_StateItem_Flag)     :: itemFlag
  integer                       :: n, i, j, i1, j1
  integer                       :: lbnd1,lbnd2
  logical                       :: do_sum_loc
  real(ESMF_KIND_R8), pointer   :: dataPtr1d(:)
  real(ESMF_KIND_R8), pointer   :: dataPtr2d(:,:)
  character(len=*)  , parameter :: subname='(MOM_cap_methods:state_getimport_2d)'
  ! ----------------------------------------------

  rc = ESMF_SUCCESS
  if (present(do_sum)) then
     do_sum_loc = do_sum
  else
     do_sum_loc = .false.
  endif

  call ESMF_StateGet(State, trim(fldname), itemFlag, rc=rc)
  if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then

    if (geomtype == ESMF_GEOMTYPE_MESH) then

      ! get field pointer
      if (present(esmf_ind)) then
         call state_getfldptr(state, trim(fldname), dataptr2d, rc)
         dataptr1d => dataptr2d(esmf_ind,:)
      else
         call state_getfldptr(state, trim(fldname), dataptr1d, rc)
      endif
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! determine output array and apply area correction if present
      n = 0
      do j = jsc,jec
        do i = isc,iec
          n = n + 1
          if (do_sum_loc) then
            if (present(areacor)) then
              output(i,j)  = output(i,j) + dataPtr1d(n) * areacor(n)
            else
              output(i,j)  = output(i,j) + dataPtr1d(n)
            endif
          else
            if (present(areacor)) then
              output(i,j)  = dataPtr1d(n) * areacor(n)
            else
              output(i,j)  = dataPtr1d(n)
            endif
          endif
        enddo
      enddo

    elseif (geomtype == ESMF_GEOMTYPE_GRID) then

      call state_getfldptr(state, trim(fldname), dataptr2d, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      lbnd1 = lbound(dataPtr2d,1)
      lbnd2 = lbound(dataPtr2d,2)

      do j = jsc, jec
        j1 = j + lbnd2 - jsc
        do i = isc, iec
          i1 = i + lbnd1 - isc
          if (do_sum_loc) then
            output(i,j) = output(i,j) + dataPtr2d(i1,j1)
          else
            output(i,j) = dataPtr2d(i1,j1)
          endif
        enddo
      enddo

    endif

  endif

end subroutine State_GetImport_2d

!> Map 3d import state field to output array (where 3rd dim is an ungridded dimension)
subroutine State_GetImport_3d(state, fldname, isc, iec, jsc, jec, lbd, ubd, output, do_sum, areacor, rc)
  type(ESMF_State)    , intent(in)    :: state   !< ESMF state
  character(len=*)    , intent(in)    :: fldname !< Field name
  integer             , intent(in)    :: isc     !< The start i-index of cell centers within
                                                 !! the computational domain
  integer             , intent(in)    :: iec     !< The end i-index of cell centers within the
                                                 !! computational domain
  integer             , intent(in)    :: jsc     !< The start j-index of cell centers within
                                                 !! the computational domain
  integer             , intent(in)    :: jec     !< The end j-index of cell centers within
                                                 !! the computational domain
  integer             , intent(in)    :: lbd     !< lower bound of ungridded dimension
  integer             , intent(in)    :: ubd     !< upper bound of ungridded dimension
  real (ESMF_KIND_R8) , intent(inout) :: output(isc:iec,jsc:jec,lbd:ubd)!< Output 3D array
  logical, optional   , intent(in)    :: do_sum  !< If true, sums the data
  real (ESMF_KIND_R8), optional,  intent(in) :: areacor(:) !< flux area correction factors
                                                           !! applicable to meshes
  integer             , intent(out)   :: rc      !< Return code

  ! local variables
  type(ESMF_StateItem_Flag)     :: itemFlag
  integer                       :: n, i, j, i1, j1, u
  integer                       :: lbnd1,lbnd2
  logical                       :: do_sum_loc
  real(ESMF_KIND_R8), pointer   :: dataPtr2d(:,:)
  character(len=*)  , parameter :: subname='(MOM_cap_methods:state_getimport_3d)'
  ! ----------------------------------------------

  rc = ESMF_SUCCESS
  if (present(do_sum)) then
     do_sum_loc = do_sum
  else
     do_sum_loc = .false.
  endif

  call ESMF_StateGet(State, trim(fldname), itemFlag, rc=rc)
  if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then

     if (geomtype == ESMF_GEOMTYPE_MESH) then

        ! get field pointer
        call state_getfldptr(state, trim(fldname), dataptr2d, rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        ! determine output array and apply area correction if present
        do u = lbd, ubd ! ungridded dims
           n = 0
           do j = jsc,jec
              do i = isc,iec
                 n = n + 1
                 if (do_sum_loc) then
                    if (present(areacor)) then
                       output(i,j,u)  = output(i,j,u) + dataPtr2d(u,n) * areacor(n)
                    else
                       output(i,j,u)  = output(i,j,u) + dataPtr2d(u,n)
                    endif
                 else
                    if (present(areacor)) then
                       output(i,j,u)  = dataPtr2d(u,n) * areacor(n)
                    else
                       output(i,j,u)  = dataPtr2d(u,n)
                    endif
                 endif
              enddo
           enddo
         enddo

     else if (geomtype == ESMF_GEOMTYPE_GRID) then
        call ESMF_LogWrite(trim(subname)//": ERROR ungridded dimensions not supported in MOM6 nuopc cap when "// &
         "ESMF_GEOMTYPE_GRID is used. Use ESMF_GEOMTYPE_MESH instead.", ESMF_LOGMSG_ERROR)
        rc = ESMF_FAILURE
        return
     endif

  endif

end subroutine State_GetImport_3d

!> Map input array to export state
subroutine State_SetExport(state, fldname, isc, iec, jsc, jec, input, ocean_grid, areacor, rc)
  type(ESMF_State)      , intent(inout) :: state   !< ESMF state
  character(len=*)      , intent(in)    :: fldname !< Field name
  integer             , intent(in)      :: isc     !< The start i-index of cell centers within
                                                   !! the computational domain
  integer             , intent(in)      :: iec     !< The end i-index of cell centers within the
                                                   !! computational domain
  integer             , intent(in)      :: jsc     !< The start j-index of cell centers within
                                                   !! the computational domain
  integer             , intent(in)      :: jec     !< The end j-index of cell centers within
                                                   !! the computational domain
  real (ESMF_KIND_R8)   , intent(in)    :: input(isc:iec,jsc:jec)!< Input 2D array
  type(ocean_grid_type) , intent(in)    :: ocean_grid !< Ocean horizontal grid
  real (ESMF_KIND_R8), optional,  intent(in) :: areacor(:) !< flux area correction factors
                                                           !! applicable to meshes
  integer               , intent(out)   :: rc         !< Return code

  ! local variables
  type(ESMF_StateItem_Flag)     :: itemFlag
  integer                       :: n, i, j, k, i1, j1, ig,jg
  integer                       :: lbnd1,lbnd2
  real(ESMF_KIND_R8), pointer   :: dataPtr1d(:)
  real(ESMF_KIND_R8), pointer   :: dataPtr2d(:,:)
  character(len=*)  , parameter :: subname='(MOM_cap_methods:state_setexport)'
  ! ----------------------------------------------

  rc = ESMF_SUCCESS

  ! Indexing notes:
  ! input array from "ocean_public" uses local indexing without halos
  ! mask from "ocean_grid" uses local indexing with halos

  call ESMF_StateGet(State, trim(fldname), itemFlag, rc=rc)
  if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then

    if (geomtype == ESMF_GEOMTYPE_MESH) then

      call state_getfldptr(state, trim(fldname), dataptr1d, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      n = 0
      do j = jsc, jec
        jg = j + ocean_grid%jsc - jsc
        do i = isc, iec
          ig = i + ocean_grid%isc - isc
          n = n+1
          dataPtr1d(n) = input(i,j) * ocean_grid%mask2dT(ig,jg)
        enddo
      enddo
      if (present(areacor)) then
        do n = 1,(size(dataPtr1d))
          dataPtr1d(n) = dataPtr1d(n) * areacor(n)
        enddo
      endif

      ! if a maskmap is provided, set exports of all eliminated cells to zero.
      if (associated(ocean_grid%Domain%maskmap)) then
        do k = n+1, size(dataPtr1d)
          dataPtr1d(k) = 0.0
        enddo
      endif

    else if (geomtype == ESMF_GEOMTYPE_GRID) then

      call state_getfldptr(state, trim(fldname), dataptr2d, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      lbnd1 = lbound(dataPtr2d,1)
      lbnd2 = lbound(dataPtr2d,2)

      do j = jsc, jec
        j1 = j + lbnd2 - jsc
        jg = j + ocean_grid%jsc - jsc
        do i = isc, iec
          i1 = i + lbnd1 - isc
          ig = i + ocean_grid%isc - isc
          dataPtr2d(i1,j1)  = input(i,j) * ocean_grid%mask2dT(ig,jg)
        enddo
      enddo

    endif

  endif

end subroutine State_SetExport

!> This subroutine writes the minimum, maximum and sum of each field
!! contained within an ESMF state.
subroutine state_diagnose(State, string, rc)

  ! ----------------------------------------------
  ! Diagnose status of State
  ! ----------------------------------------------

  type(ESMF_State), intent(in)    :: state  !< An ESMF State
  character(len=*), intent(in)    :: string !< A string indicating whether the State is an
                                            !! import or export State
  integer         , intent(out)   :: rc     !< Return code

  ! local variables
  integer                         :: i,j,n
  type(ESMf_Field)                :: lfield
  integer                         :: fieldCount, lrank
  character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
  real(ESMF_KIND_R8), pointer     :: dataPtr1d(:)
  real(ESMF_KIND_R8), pointer     :: dataPtr2d(:,:)
  character(len=*),parameter      :: subname='(state_diagnose)'
  character(len=ESMF_MAXSTR)      :: msgString
  ! ----------------------------------------------

  call ESMF_StateGet(state, itemCount=fieldCount, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  allocate(lfieldnamelist(fieldCount))

  call ESMF_StateGet(state, itemNameList=lfieldnamelist, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  do n = 1, fieldCount

    call ESMF_StateGet(state, itemName=lfieldnamelist(n), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call field_getfldptr(lfield, fldptr1=dataPtr1d, fldptr2=dataPtr2d, rank=lrank, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (lrank == 0) then
      ! no local data
    elseif (lrank == 1) then
      if (size(dataPtr1d) > 0) then
        write(msgString,'(A,3g14.7,i8)') trim(string)//': '//trim(lfieldnamelist(n)), &
              minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
      else
        write(msgString,'(A,a)') trim(string)//': '//trim(lfieldnamelist(n))," no data"
      endif
    elseif (lrank == 2) then
      if (size(dataPtr2d) > 0) then
        write(msgString,'(A,3g14.7,i8)') trim(string)//': '//trim(lfieldnamelist(n)), &
              minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
      else
        write(msgString,'(A,a)') trim(string)//': '//trim(lfieldnamelist(n))," no data"
      endif
    else
      call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
    endif
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
  enddo

  deallocate(lfieldnamelist)

end subroutine state_diagnose

!> Obtain a pointer to a rank 1 or rank 2 ESMF field
subroutine field_getfldptr(field, fldptr1, fldptr2, rank, abort, rc)

  ! input/output variables
  type(ESMF_Field)  , intent(in)                        :: field        !< An ESMF field
  real(ESMF_KIND_R8), pointer , intent(inout), optional :: fldptr1(:)   !< A pointer to a rank 1 ESMF field
  real(ESMF_KIND_R8), pointer , intent(inout), optional :: fldptr2(:,:) !< A pointer to a rank 2 ESMF field
  integer           , intent(out)            , optional :: rank         !< Field rank
  logical           , intent(in)             , optional :: abort        !< Abort code
  integer           , intent(out)            , optional :: rc           !< Return code

  ! local variables
  type(ESMF_GeomType_Flag)    :: geomtype
  type(ESMF_FieldStatus_Flag) :: status
  type(ESMF_Mesh)             :: lmesh
  integer                     :: lrank, nnodes, nelements
  logical                     :: labort
  character(len=*), parameter :: subname='(field_getfldptr)'
  ! ----------------------------------------------

  if (.not.present(rc)) then
    call ESMF_LogWrite(trim(subname)//": ERROR rc not present ", &
         ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
    rc = ESMF_FAILURE
    return
  endif

  rc = ESMF_SUCCESS

  labort = .true.
  if (present(abort)) then
    labort = abort
  endif
  lrank = -99

  call ESMF_FieldGet(field, status=status, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  if (status /= ESMF_FIELDSTATUS_COMPLETE) then
    lrank = 0
    if (labort) then
      call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE
      return
    else
      call ESMF_LogWrite(trim(subname)//": WARNING data not allocated ", ESMF_LOGMSG_INFO)
    endif
  else

    call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (geomtype == ESMF_GEOMTYPE_GRID) then
      call ESMF_FieldGet(field, rank=lrank, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    elseif (geomtype == ESMF_GEOMTYPE_MESH) then
      call ESMF_FieldGet(field, rank=lrank, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldGet(field, mesh=lmesh, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (nnodes == 0 .and. nelements == 0) lrank = 0
    else
      call ESMF_LogWrite(trim(subname)//": ERROR geomtype not supported ", &
           ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE
      return
    endif ! geomtype

    if (lrank == 0) then
      call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", &
           ESMF_LOGMSG_INFO)
    elseif (lrank == 1) then
      if (.not.present(fldptr1)) then
        call ESMF_LogWrite(trim(subname)//": ERROR missing rank=1 array ", &
             ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
        rc = ESMF_FAILURE
        return
      endif
      call ESMF_FieldGet(field, farrayPtr=fldptr1, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    elseif (lrank == 2) then
      if (.not.present(fldptr2)) then
        call ESMF_LogWrite(trim(subname)//": ERROR missing rank=2 array ", &
             ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
        rc = ESMF_FAILURE
        return
      endif
      call ESMF_FieldGet(field, farrayPtr=fldptr2, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//": ERROR in rank ", &
           ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
      rc = ESMF_FAILURE
      return
    endif

  endif  ! status

  if (present(rank)) then
    rank = lrank
  endif

end subroutine field_getfldptr

!> Returns true if ESMF_LogFoundError() determines that rc is an error code. Otherwise false.
logical function ChkErr(rc, line, file)
  integer, intent(in) :: rc            !< return code to check
  integer, intent(in) :: line          !< Integer source line number
  character(len=*), intent(in) :: file !< User-provided source file name
  integer :: lrc
  ChkErr = .false.
  lrc = rc
  if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
    ChkErr = .true.
  endif
end function ChkErr

end module MOM_cap_methods
