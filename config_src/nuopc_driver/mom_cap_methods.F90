!> This is the main driver for MOM6 in CIME
module mom_cap_methods

! This file is part of MOM6. See LICENSE.md for the license.

! mct modules
use ESMF
use perf_mod,        only: t_startf, t_stopf
use ocean_model_mod, only: ocean_public_type, ocean_state_type
use ocean_model_mod, only: ice_ocean_boundary_type
use MOM_grid,        only: ocean_grid_type
use MOM_domains,     only: pass_var
use mpp_domains_mod, only: mpp_get_compute_domain

! By default make data private
implicit none; private

! Public member functions
public :: ocn_export
public :: ocn_import

integer :: rc,dbrc
character(len=1024) :: tmpstr

!---------------------------
contains
!---------------------------

!> Maps outgoing ocean data to ESMF State
!! See \ref section_ocn_export for a summary of the data
!! that is transferred from MOM6 to MCT.
subroutine ocn_export(ocean_public, grid, exportState)
  type(ocean_public_type), intent(in)    :: ocean_public !< Ocean surface state
  type(ocean_grid_type),   intent(in)    :: grid       !< Ocean model grid
  type(ESMF_State),        intent(inout) :: exportState !< outgoing data
  ! Local variables
  real, dimension(grid%isd:grid%ied,grid%jsd:grid%jed) :: ssh !< Local copy of sea_lev with updated halo
  integer :: i, j, i1, j1, isc, iec, jsc, jec  !< Grid indices
  integer :: lbnd1, lbnd2, ubnd1, ubnd2
  real :: slp_L, slp_R, slp_C, slope, u_min, u_max
  real(ESMF_KIND_R8), pointer :: dataPtr_omask(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_t(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_s(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_u(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_v(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_q(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_dhdx(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_dhdy(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_bldepth(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_fswpen(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_roce_16O(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_roce_HDO(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_fco2_ocn(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_fdms_ocn(:,:)

  character(len=*),parameter :: subname = '(ocn_export)'

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
  call State_getFldPtr(exportState,"Fioo_q", dataPtr_q, rc=rc)
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
  call State_getFldPtr(exportState,"So_bldepth", dataPtr_bldepth, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call State_getFldPtr(exportState,"So_fswpen", dataPtr_fswpen, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
!  call State_getFldPtr(exportState,"So_roce_16O", dataPtr_roce_16O, rc=rc)
!  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!    line=__LINE__, &
!    file=__FILE__)) &
!    return  ! bail out
!  call State_getFldPtr(exportState,"So_roce_HDO", dataPtr_roce_HDO, rc=rc)
!  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!    line=__LINE__, &
!    file=__FILE__)) &
!    return  ! bail out
!  call State_getFldPtr(exportState,"Faoo_fco2_ocn", dataPtr_fco2_ocn, rc=rc)
!  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!    line=__LINE__, &
!    file=__FILE__)) &
!    return  ! bail out
!  call State_getFldPtr(exportState,"Faoo_fdms_ocn", dataPtr_fdms_ocn, rc=rc)
!  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!    line=__LINE__, &
!    file=__FILE__)) &
!    return  ! bail out

  lbnd1 = lbound(dataPtr_t,1)
  ubnd1 = ubound(dataPtr_t,1)
  lbnd2 = lbound(dataPtr_t,2)
  ubnd2 = ubound(dataPtr_t,2)

  call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

!tcx
  write(tmpstr,'(a,6i8)') subname//'tcx1',lbnd1,ubnd1,lbnd2,ubnd2
  call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
  write(tmpstr,'(a,6i8)') subname//'tcx2',isc,iec,jsc,jec
  call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
  write(tmpstr,'(a,6i8)') subname//'tcx3',lbound(ssh,1),ubound(ssh,1),lbound(ssh,2),ubound(ssh,2)
  call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
  write(tmpstr,'(a,6i8)') subname//'tcx4',lbound(ocean_public%sea_lev,1),ubound(ocean_public%sea_lev,1),lbound(ocean_public%sea_lev,2),ubound(ocean_public%sea_lev,2)
  call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
  write(tmpstr,'(a,6i8)') subname//'tcx5',lbound(grid%mask2dT,1),ubound(grid%mask2dT,1),lbound(grid%mask2dT,2),ubound(grid%mask2dT,2)
  call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

  ! Copy from ocean_public to exportstate. ocean_public uses global indexing with no halos.
  ! The mask comes from "grid" that uses the usual MOM domain that has halos
  ! and does not use global indexing.
  do j = jsc, jec
    j1 = j + lbnd2 - jsc
    do i = isc, iec
      i1 = i + lbnd1 - isc
      ! surface temperature in Kelvin
      dataPtr_t(i1,j1) = ocean_public%t_surf(i,j)  * grid%mask2dT(i,j)
      dataPtr_s(i1,j1) = ocean_public%s_surf(i,j) * grid%mask2dT(i,j)
      dataPtr_u(i1,j1) = (grid%cos_rot(i,j) * ocean_public%u_surf(i,j) &
                        - grid%sin_rot(i,j) * ocean_public%v_surf(i,j)) * grid%mask2dT(i,j)
      dataPtr_v(i1,j1) = (grid%cos_rot(i,j) * ocean_public%v_surf(i,j) &
                        + grid%sin_rot(i,j) * ocean_public%u_surf(i,j)) * grid%mask2dT(i,j)
      ! Make a copy of ssh in order to do a halo update. We use the usual MOM domain
      ! in order to update halos. i.e. does not use global indexing.
      ssh(i,j) = ocean_public%sea_lev(i,j)
    end do
  end do

  ! Update halo of ssh so we can calculate gradients
  call pass_var(ssh, grid%domain)

  ! d/dx ssh
  do j=jsc, jec
    j1 = j + lbnd2 - jsc
    do i=isc,iec
      i1 = i + lbnd1 - isc
      ! This is a simple second-order difference
      !dataPtr_dhdx(i1,j1) = 0.5 * (ssh(i+1,j) - ssh(i-1,j)) * grid%IdxT(i,j) * grid%mask2dT(i,j)
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
      end if
      dataPtr_dhdx(i1,j1) = slope * grid%IdxT(i,j) * grid%mask2dT(i,j)
      if (grid%mask2dT(i,j)==0.) dataPtr_dhdx(i1,j1) = 0.0
    end do
  end do

  ! d/dy ssh
  do j=jsc, jec
    j1 = j + lbnd2 - jsc
    do i=isc,iec
      i1 = i + lbnd1 - isc
      ! This is a simple second-order difference
      !dataPtr_dhdy(i1,j1) = 0.5 * (ssh(i,j+1) - ssh(i,j-1)) * grid%IdyT(i,j) * grid%mask2dT(i,j)
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
      end if
      dataPtr_dhdy(i1,j1) = slope * grid%IdyT(i,j) * grid%mask2dT(i,j)
      if (grid%mask2dT(i,j)==0.) dataPtr_dhdy(i1,j1) = 0.0
    end do
  end do

end subroutine ocn_export


!> This function has a few purposes: 1) it allocates and initializes the data
!! in the fluxes structure; 2) it imports surface fluxes using data from
!! the coupler; and 3) it can apply restoring in SST and SSS.
!! See \ref section_ocn_import for a summary of the surface fluxes that are
!! passed from MCT to MOM6, including fluxes that need to be included in
!! the future.
subroutine ocn_import(ocean_public, grid, importState, ice_ocean_boundary)
  type(ocean_public_type), intent(in)    :: ocean_public !< Ocean surface state
  type(ocean_grid_type),   intent(in)    :: grid       !< Ocean model grid
  type(ESMF_State),        intent(inout) :: importState !< incoming data
  type(ice_ocean_boundary_type), intent(inout) :: ice_ocean_boundary !< Ocean boundary forcing

  integer :: i, j, i1, j1, isc, iec, jsc, jec  !< Grid indices
  real(ESMF_KIND_R8) :: c1,c2,c3,c4
  integer :: lbnd1, lbnd2, ubnd1, ubnd2

  real(ESMF_KIND_R8), pointer :: dataPtr_p(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_ifrac(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_duu10n(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_co2prog(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_co2diag(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_lamult(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_ustokes(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_vstokes(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_taux(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_tauy(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_sen(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_lat(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_evap(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_osalt(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_lwdn(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_lwup(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_swnet(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_rofl(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_rofi(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_meltw(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_melth(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_iosalt(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_prec(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_rain(:,:)
  real(ESMF_KIND_R8), pointer :: dataPtr_snow(:,:)

  character(len=*),parameter :: subname = '(ocn_import)'

  call State_getFldPtr(importState,'Sa_pslv', dataPtr_p,rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call State_getFldPtr(importState,'Si_ifrac', dataPtr_ifrac,rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call State_getFldPtr(importState,"So_duu10n", dataPtr_duu10n, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
!  call State_getFldPtr(importState,"Sa_co2prog", dataPtr_co2prog, rc=rc)
!  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!    line=__LINE__, &
!    file=__FILE__)) &
!    return  ! bail out
!  call State_getFldPtr(importState,"Sa_co2diag", dataPtr_co2diag, rc=rc)
!  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!    line=__LINE__, &
!    file=__FILE__)) &
!    return  ! bail out
  call State_getFldPtr(importState,"Sw_lamult", dataPtr_lamult, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call State_getFldPtr(importState,"Sw_ustokes", dataPtr_ustokes, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call State_getFldPtr(importState,"Sw_vstokes", dataPtr_vstokes, rc=rc)
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
  call State_getFldPtr(importState,"Foxx_salt" , dataPtr_osalt, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call State_getFldPtr(importState,"Foxx_lwdn" , dataPtr_lwdn, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call State_getFldPtr(importState,"Foxx_lwup" , dataPtr_lwup, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call State_getFldPtr(importState,"Foxx_swnet", dataPtr_swnet, rc=rc)
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
  call State_getFldPtr(importState,"Foxx_meltw", dataPtr_meltw, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call State_getFldPtr(importState,"Foxx_melth", dataPtr_melth, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call State_getFldPtr(importState,"Foxx_salt" , dataPtr_iosalt, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call State_getFldPtr(importState,"Foxx_prec" , dataPtr_prec, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call State_getFldPtr(importState,"Foxx_rain" , dataPtr_rain, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call State_getFldPtr(importState,"Foxx_snow" , dataPtr_snow, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out

  lbnd1 = lbound(dataPtr_p,1)
  ubnd1 = ubound(dataPtr_p,1)
  lbnd2 = lbound(dataPtr_p,2)
  ubnd2 = ubound(dataPtr_p,2)

  call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

!tcx
!  write(tmpstr,'(a,6i8)') subname//'tcx1',lbnd1,ubnd1,lbnd2,ubnd2
!  call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
!  write(tmpstr,'(a,6i8)') subname//'tcx2',isc,iec,jsc,jec
!  call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
!  write(tmpstr,'(a,6i8)') subname//'tcx3',i,j,i1,j1
!  call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
!  write(tmpstr,'(a,6i8)') subname//'tcx4',lbound(ice_ocean_boundary%p,1),ubound(ice_ocean_boundary%p,1),lbound(ice_ocean_boundary%p,2),ubound(ice_ocean_boundary%p,2)
!  call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
!  write(tmpstr,'(a,6i8)') subname//'tcx5',lbound(grid%mask2dT,1),ubound(grid%mask2dT,1),lbound(grid%mask2dT,2),ubound(grid%mask2dT,2)
!  call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

  do j = jsc,jec
    do i = isc,iec
      i1 = i + lbnd1 - isc
      j1 = j + lbnd2 - jsc

      ice_ocean_boundary%p(i,j) = GRID%mask2dT(i,j) * dataPtr_p(i1,j1)

      ice_ocean_boundary%u_flux(i,j) = (GRID%cos_rot(i1,j1)*dataPtr_taux(i1,j1) + &
                                        GRID%sin_rot(i1,j1)*dataPtr_tauy(i1,j1))
      ice_ocean_boundary%v_flux(i,j) = (GRID%cos_rot(i1,j1)*dataPtr_tauy(i1,j1) + &
                                        GRID%sin_rot(i1,j1)*dataPtr_taux(i1,j1))

      ice_ocean_boundary%t_flux(i,j)  = dataPtr_sen(i1,j1) * GRID%mask2dT(i,j)
      ice_ocean_boundary%q_flux(i,j)  = dataPtr_evap(i1,j1) * GRID%mask2dT(i,j)
!     ice_ocean_boundary%latent(i,j)  = dataPtr_lat(i1,j1) * GRID%mask2dT(i,j)
      ice_ocean_boundary%lw_flux(i,j) = (dataPtr_lwup(i1,j1) + dataPtr_lwdn(i1,j1)) * GRID%mask2dT(i,j)

! tcx TO DO c1-c4
      c1 = 0.25_ESMF_KIND_R8
      c2 = 0.25_ESMF_KIND_R8
      c3 = 0.25_ESMF_KIND_R8
      c4 = 0.25_ESMF_KIND_R8
      ice_ocean_boundary%sw_flux_vis_dir(i,j) = GRID%mask2dT(i,j) * dataPtr_swnet(i1,j1)*c1
      ice_ocean_boundary%sw_flux_vis_dif(i,j) = GRID%mask2dT(i,j) * dataPtr_swnet(i1,j1)*c2
      ice_ocean_boundary%sw_flux_nir_dir(i,j) = GRID%mask2dT(i,j) * dataPtr_swnet(i1,j1)*c3
      ice_ocean_boundary%sw_flux_nir_dif(i,j) = GRID%mask2dT(i,j) * dataPtr_swnet(i1,j1)*c4

!      ice_ocean_boundary%sw(i,j) = ice_ocean_boundary%sw_flux_vis_dir(i,j) + ice_ocean_boundary%sw_flux_vis_dif(i,j) + &
!                       ice_ocean_boundary%sw_flux_nir_dir(i,j) + ice_ocean_boundary%sw_flux_nir_dif(i,j)

      ice_ocean_boundary%lprec(i,j) = dataPtr_rain(i1,j1) * GRID%mask2dT(i,j)
      ice_ocean_boundary%fprec(i,j) = dataPtr_snow(i1,j1) * GRID%mask2dT(i,j)
      ice_ocean_boundary%runoff(i,j) = (dataPtr_rofl(i1,j1)+dataPtr_rofi(i1,j1)) * GRID%mask2dT(i,j)
      ice_ocean_boundary%runoff_hflx(i,j)  = 0.0 * GRID%mask2dT(i,j)
      ice_ocean_boundary%calving(i,j)      = 0.0 * GRID%mask2dT(i,j)
      ice_ocean_boundary%calving_hflx(i,j) = 0.0 * GRID%mask2dT(i,j)
      ice_ocean_boundary%salt_flux(i,j) = GRID%mask2dT(i,j)*dataPtr_iosalt(i1,j1)
      ice_ocean_boundary%salt_flux(i,j) = GRID%mask2dT(i,j)*(dataPtr_osalt(i1,j1) + ice_ocean_boundary%salt_flux(i,j))

      ice_ocean_boundary%ustar_berg(i,j) = 0.0 * GRID%mask2dT(i,j)
      ice_ocean_boundary%area_berg(i,j)  = 0.0 * GRID%mask2dT(i,j)
      ice_ocean_boundary%mass_berg(i,j)  = 0.0 * GRID%mask2dT(i,j)
      ice_ocean_boundary%mi(i,j)         = 0.0 * GRID%mask2dT(i,j)

    enddo
  enddo

end subroutine ocn_import


  !-----------------------------------------------------------------------------                                                                      

  subroutine State_GetFldPtr(ST, fldname, fldptr, rc)
    type(ESMF_State), intent(in) :: ST
    character(len=*), intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:,:)
    integer, intent(out), optional :: rc

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
