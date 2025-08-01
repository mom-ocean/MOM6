!> Provides diagnostics of work due to a given diffusivity
module MOM_diagnose_kdwork

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : diag_ctrl, time_type, post_data, register_diag_field
use MOM_diag_mediator, only : register_scalar_field
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_spatial_means, only : global_area_integral

implicit none ; private

#include <MOM_memory.h>

public vbf_CS
public kdwork_diagnostics
public Allocate_VBF_CS
public Deallocate_VBF_CS
public KdWork_init
public KdWork_end

!> This structure has memory for used in calculating diagnostics of diffusivity
!! many of the diffusivity diagnostics are copies of other 3d arrays.  It could
!! be written more efficiently, but it is less intrusive to copy into this structure
!! and do all calculations in this module.  These diagnostics may be expensive for
!! routine use.
type vbf_CS
  ! 3d varying Kd contributions
  real, pointer, dimension(:,:,:) :: &
    Bflx_salt => NULL(), & !< Salinity contribution to buoyancy flux at interfaces
                           !! [H Z T-3 ~> m2 s-3 or kg m-1 s-3 = W m-3]
    Bflx_temp => NULL(), & !< Temperature contribution to buoyancy flux at interfaces
                           !! [H Z T-3 ~> m2 s-3 or kg m-1 s-3 = W m-3]
    Bflx_salt_dz => NULL(), & !< Salinity contribution to integral of buoyancy flux over layer
                              !! [H Z2 T-3 ~> m3 s-3 or kg m-1 s-3 = W m-2]
    Bflx_temp_dz => NULL(), & !< Temperature contribution to integral of buoyancy flux over layer
                              !! [H Z2 T-3 ~> m3 s-3 or kg m-1 s-3 = W m-2]
    ! The following are all allocatable arrays that store copies of process driven Kd, so that
    ! the process driven buoyancy flux and work can be derived at the end of the time step.
    Kd_salt => NULL(), &   !< total diapycnal diffusivity of salt at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_temp => NULL(), &   !< total diapycnal diffusivity of heat at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_BBL => NULL(), &    !< diapycnal diffusivity due to BBL at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_ePBL => NULL(), &   !< diapycnal diffusivity due to ePBL at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_KS => NULL(), &     !< diapycnal diffusivity due to Kappa Shear at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_bkgnd => NULL(), &  !< diapycnal diffusivity due to Kd_bkgnd at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_ddiff_S => NULL(), &!< diapycnal diffusivity due to double diffusion of salt at interfaces
                           !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_ddiff_T => NULL(), &!< diapycnal diffusivity due to double diffusion of heat at interfaces
                           !![H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_leak => NULL(), &   !< diapycnal diffusivity due to Kd_leak at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_quad => NULL(), &   !< diapycnal diffusivity due to Kd_quad at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_itidal => NULL(), & !< diapycnal diffusivity due to Kd_itidal at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_Froude => NULL(), & !< diapycnal diffusivity due to Kd_Froude at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_slope => NULL(), &  !< diapycnal diffusivity due to Kd_slope at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_lowmode => NULL(), &!< diapycnal diffusivity due to Kd_lowmode at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_Niku => NULL(), &   !< diapycnal diffusivity due to Kd_Niku at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_itides => NULL()    !< diapycnal diffusivity due to Kd_itides at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]

  ! Constant Kd contributions
  real :: Kd_add !< spatially uniform additional diapycnal diffusivity at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
                 !! a diagnostic for this diffusivity is not yet included, but this makes it straightforward to add

  !>@{ Diagnostic IDs
  integer :: id_Bdif     = -1, id_Bdif_salt  = -1, id_Bdif_temp  = -1
  integer :: id_Bdif_dz  = -1, id_Bdif_salt_dz  = -1, id_Bdif_temp_dz  = -1
  integer :: id_Bdif_idz  = -1, id_Bdif_salt_idz  = -1, id_Bdif_temp_idz  = -1
  integer :: id_Bdif_idV  = -1, id_Bdif_salt_idV  = -1, id_Bdif_temp_idV  = -1
  integer :: id_Bdif_ePBL  = -1, id_Bdif_dz_ePBL  = -1, id_Bdif_idz_ePBL  = -1, id_Bdif_idV_ePBL  = -1
  integer :: id_Bdif_BBL  = -1, id_Bdif_dz_BBL  = -1, id_Bdif_idz_BBL  = -1, id_Bdif_idV_BBL  = -1
  integer :: id_Bdif_KS  = -1, id_Bdif_dz_KS  = -1, id_Bdif_idz_KS  = -1, id_Bdif_idV_KS  = -1
  integer :: id_Bdif_bkgnd  = -1, id_Bdif_dz_bkgnd  = -1, id_Bdif_idz_bkgnd  = -1, id_Bdif_idV_bkgnd  = -1
  integer :: id_Bdif_ddiff_temp  = -1, id_Bdif_ddiff_salt  = -1
  integer :: id_Bdif_dz_ddiff_temp  = -1, id_Bdif_dz_ddiff_salt  = -1
  integer :: id_Bdif_idz_ddiff_temp  = -1, id_Bdif_idz_ddiff_salt  = -1
  integer :: id_Bdif_idV_ddiff_temp  = -1, id_Bdif_idV_ddiff_salt  = -1
  integer :: id_Bdif_leak  = -1, id_Bdif_dz_leak  = -1, id_Bdif_idz_leak  = -1, id_Bdif_idV_leak  = -1
  integer :: id_Bdif_quad  = -1, id_Bdif_dz_quad  = -1, id_Bdif_idz_quad  = -1, id_Bdif_idV_quad  = -1
  integer :: id_Bdif_itidal  = -1, id_Bdif_dz_itidal  = -1, id_Bdif_idz_itidal  = -1, id_Bdif_idV_itidal  = -1
  integer :: id_Bdif_Froude  = -1, id_Bdif_dz_Froude  = -1, id_Bdif_idz_Froude  = -1, id_Bdif_idV_Froude  = -1
  integer :: id_Bdif_slope  = -1, id_Bdif_dz_slope  = -1, id_Bdif_idz_slope  = -1, id_Bdif_idV_slope  = -1
  integer :: id_Bdif_lowmode  = -1, id_Bdif_dz_lowmode  = -1, id_Bdif_idz_lowmode  = -1, id_Bdif_idV_lowmode  = -1
  integer :: id_Bdif_Niku  = -1, id_Bdif_dz_Niku  = -1, id_Bdif_idz_Niku  = -1, id_Bdif_idV_Niku  = -1
  integer :: id_Bdif_itides  = -1, id_Bdif_dz_itides  = -1, id_Bdif_idz_itides  = -1, id_Bdif_idV_itides  = -1
  !>@}

  logical :: do_bflx_salt = .false.  !< Logical flag to indicate if N2_salt should be computed
  logical :: do_bflx_temp = .false.  !< Logical flag to indicate if N2_temp should be computed
  logical :: do_bflx_salt_dz = .false.  !< Logical flag to indicate if N2_salt should be computed
  logical :: do_bflx_temp_dz = .false.  !< Logical flag to indicate if N2_temp should be computed

end type vbf_CS

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> Loop over all implemented diffusivities to diagnose and output Kd Work/buoyancy fluxes
subroutine KdWork_Diagnostics(G,GV,US,diag,VBF,N2_Salt,N2_Temp,dz)
  type(ocean_grid_type),      intent(in)    :: G       !< Grid type
  type(verticalGrid_type),    intent(in)    :: GV      !< ocean vertical grid structure
  type(unit_scale_type),      intent(in)    :: US      !< A dimensional unit scaling type
  type(diag_ctrl), target,    intent(inout) :: diag    !< regulates diagnostic output
  type (vbf_CS),              intent(inout) :: VBF     !< Vertical buoyancy flux structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                              intent(in)    :: N2_Salt !< Buoyancy frequency [T-2 ~> s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                              intent(in)    :: N2_Temp !< Buoyancy frequency [T-2 ~> s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in)    :: dz      !< Grid spacing [Z ~> m]

  ! Work arrays for computing buoyancy flux integrals
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: work3d_i
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: work3d_l
  real, dimension(SZI_(G),SZJ_(G)) :: work2d, work2d_salt, work2d_temp
  real :: work, work_salt, work_temp

  integer :: i, j, k, nz, isc, iec, jsc, jec

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ;

  nz = GV%ke

  ! Compute total fluxes
  if (VBF%id_Bdif_dz>0 .or. VBF%id_Bdif_salt_dz>0 .or. VBF%id_Bdif_temp_dz>0 .or. &
      VBF%id_Bdif_idz>0 .or. VBF%id_Bdif_salt_idz>0 .or. VBF%id_Bdif_temp_idz>0 .or. &
      VBF%id_Bdif_idV>0 .or. VBF%id_Bdif_salt_idV>0 .or. VBF%id_Bdif_temp_idV>0 ) then ! Doing vertical integrals
    ! Do Salt
    if (VBF%id_Bdif_salt_dz>0 .or. VBF%id_Bdif_dz>0 .or. VBF%id_Bdif_salt>0 .or. VBF%id_Bdif>0 .or. &
        VBF%id_Bdif_idz>0 .or. VBF%id_Bdif_salt_idz>0 .or. VBF%id_Bdif_idV>0 .or. VBF%id_Bdif_salt_idV>0) &
      call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_salt, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    ! Do Temp
    if (VBF%id_Bdif_temp_dz>0 .or. VBF%id_Bdif_dz>0 .or. VBF%id_Bdif_temp>0 .or. VBF%id_Bdif>0 .or. &
        VBF%id_Bdif_idz>0 .or. VBF%id_Bdif_temp_idz>0 .or. VBF%id_Bdif_idV>0 .or. VBF%id_Bdif_temp_idV>0) &
       call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_temp, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_temp_idz>0 .or. VBF%id_Bdif_idz>0) then
      work2d_temp(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d_temp(i,j) = work2d_temp(i,j) + VBF%Bflx_temp_dz(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_temp_idV>0 .or. VBF%id_Bdif_idV>0) then
      work_temp = 0.0
      do k = 1,nz
        work_temp = work_temp + global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, &
                    tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
      enddo
    endif
    if (VBF%id_Bdif_salt_idz>0 .or. VBF%id_Bdif_idz>0) then
      work2d_salt(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d_salt(i,j) = work2d_salt(i,j) + VBF%Bflx_salt_dz(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_salt_idV>0 .or. VBF%id_Bdif_idV>0) then
      work_salt = 0.0
      do k = 1,nz
        work_salt = work_salt + global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, &
                    tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
      enddo
    endif
    work = work_temp + work_salt
    do j = jsc,jec ; do i = isc,iec
      work2d(i,j) = work2d_temp(i,j) + work2d_salt(i,j)
    enddo ; enddo
  elseif (VBF%id_Bdif>0 .or. VBF%id_Bdif_salt>0 .or. VBF%id_Bdif_temp>0) then ! Not doing vertical integrals
    ! Do Salt
    if (VBF%id_Bdif_salt>0 .or. VBF%id_Bdif>0) &
      call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_salt, VBF%Bflx_salt)
    if (VBF%id_Bdif_temp>0 .or. VBF%id_Bdif>0) &
      call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_temp, VBF%Bflx_temp)
  endif
  ! Post total fluxes
  if (VBF%id_Bdif_salt>0) call post_data(VBF%id_Bdif_salt, VBF%Bflx_salt, diag)
  if (VBF%id_Bdif_temp>0) call post_data(VBF%id_Bdif_temp, VBF%Bflx_temp, diag)
  if (VBF%id_Bdif>0) then
    work3d_i(:,:,:) = 0.0
    do k = 1,nz+1 ; do j = jsc,jec ; do i = isc,iec
      work3d_i(i,j,k) = VBF%Bflx_temp(i,j,k) + VBF%Bflx_salt(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif, work3d_i, diag)
  endif
  if (VBF%id_Bdif_salt_dz>0) call post_data(VBF%id_Bdif_salt_dz, VBF%Bflx_salt_dz, diag)
  if (VBF%id_Bdif_temp_dz>0) call post_data(VBF%id_Bdif_temp_dz, VBF%Bflx_temp_dz, diag)
  if (VBF%id_Bdif_dz>0) then
    work3d_l(:,:,:) = 0.0
    do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
      work3d_l(i,j,k) = VBF%Bflx_temp_dz(i,j,k) + VBF%Bflx_salt_dz(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_dz, work3d_l, diag)
  endif
  if (VBF%id_Bdif_salt_idz>0) call post_data(VBF%id_Bdif_salt_idz, work2d_salt, diag)
  if (VBF%id_Bdif_temp_idz>0) call post_data(VBF%id_Bdif_temp_idz, work2d_temp, diag)
  if (VBF%id_Bdif_idz>0) call post_data(VBF%id_Bdif_idz, work2d, diag)
  if (VBF%id_Bdif_salt_idV>0) call post_data(VBF%id_Bdif_salt_idV, work_salt, diag)
  if (VBF%id_Bdif_temp_idV>0) call post_data(VBF%id_Bdif_temp_idV, work_temp, diag)
  if (VBF%id_Bdif_idV>0) call post_data(VBF%id_Bdif_idV, work, diag)

  ! Compute ePBL fluxes
  if (VBF%id_Bdif_dz_ePBL>0.or.VBF%id_Bdif_idz_ePBL>0.or.VBF%id_Bdif_idV_ePBL>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_ePBL, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_ePBL, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_idz_ePBL>0) then
      work2d(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d(i,j) = work2d(i,j) + (VBF%Bflx_salt_dz(i,j,k) + VBF%Bflx_temp_dz(i,j,k))
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_ePBL>0) then
      work = 0.0
      do k = 1,nz
        work = work + &
               (global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3) + &
                global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3))
      enddo
    endif
  elseif (VBF%id_Bdif_ePBL>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_ePBL, VBF%Bflx_salt)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_ePBL, VBF%Bflx_temp)
  endif
  ! Post ePBL fluxes
  if (VBF%id_Bdif_ePBL>0) then
    work3d_i(:,:,:) = 0.0
    do k = 1,nz+1 ; do j = jsc,jec ; do i = isc,iec
      work3d_i(i,j,k) = VBF%Bflx_temp(i,j,k) + VBF%Bflx_salt(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_ePBL, work3d_i, diag)
  endif
  if (VBF%id_Bdif_dz_ePBL>0) then
    work3d_l(:,:,:) = 0.0
    do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
      work3d_l(i,j,k) = VBF%Bflx_temp_dz(i,j,k) + VBF%Bflx_salt_dz(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_dz_ePBL, work3d_l, diag)
  endif
  if (VBF%id_Bdif_idz_ePBL>0) call post_data(VBF%id_Bdif_idz_ePBL, work2d, diag)
  if (VBF%id_Bdif_idV_ePBL>0) call post_data(VBF%id_Bdif_idV_ePBL, work, diag)

  ! Compute BBL fluxes
  if (VBF%id_Bdif_dz_BBL>0.or.VBF%id_Bdif_idz_BBL>0.or.VBF%id_Bdif_idV_BBL>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_BBL, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_BBL, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_idz_BBL>0) then
      work2d(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d(i,j) = work2d(i,j) + (VBF%Bflx_salt_dz(i,j,k) + VBF%Bflx_temp_dz(i,j,k))
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_BBL>0) then
      work = 0.0
      do k = 1,nz
        work = work + &
               (global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3) + &
                global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3))
      enddo
    endif
  elseif (VBF%id_Bdif_BBL>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_BBL, VBF%Bflx_salt)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_BBL, VBF%Bflx_temp)
  endif
  ! Post BBL fluxes
  if (VBF%id_Bdif_BBL>0) then
    work3d_i(:,:,:) = 0.0
    do k = 1,nz+1 ; do j = jsc,jec ; do i = isc,iec
      work3d_i(i,j,k) = VBF%Bflx_temp(i,j,k) + VBF%Bflx_salt(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_BBL, work3d_i, diag)
  endif
  if (VBF%id_Bdif_dz_BBL>0) then
    work3d_l(:,:,:) = 0.0
    do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
      work3d_l(i,j,k) = VBF%Bflx_temp_dz(i,j,k) + VBF%Bflx_salt_dz(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_dz_BBL, work3d_l, diag)
  endif
  if (VBF%id_Bdif_idz_BBL>0) call post_data(VBF%id_Bdif_idz_BBL, work2d, diag)
  if (VBF%id_Bdif_idV_BBL>0) call post_data(VBF%id_Bdif_idV_BBL, work, diag)

  ! Compute Kappa Shear fluxes
  if (VBF%id_Bdif_dz_KS>0.or.VBF%id_Bdif_idz_KS>0.or.VBF%id_Bdif_idV_KS>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_KS, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_KS, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_idz_KS>0) then
      work2d(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d(i,j) = work2d(i,j) + (VBF%Bflx_salt_dz(i,j,k) + VBF%Bflx_temp_dz(i,j,k))
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_KS>0) then
      work = 0.0
      do k = 1,nz
        work = work + &
               (global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3) + &
                global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3))
      enddo
    endif
  elseif (VBF%id_Bdif_KS>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_KS, VBF%Bflx_salt)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_KS, VBF%Bflx_temp)
  endif
  ! Post Kappa Shear fluxes
  if (VBF%id_Bdif_KS>0) then
    work3d_i(:,:,:) = 0.0
    do k = 1,nz+1 ; do j = jsc,jec ; do i = isc,iec
      work3d_i(i,j,k) = VBF%Bflx_temp(i,j,k) + VBF%Bflx_salt(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_KS, work3d_i, diag)
  endif
  if (VBF%id_Bdif_dz_KS>0) then
    work3d_l(:,:,:) = 0.0
    do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
      work3d_l(i,j,k) = VBF%Bflx_temp_dz(i,j,k) + VBF%Bflx_salt_dz(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_dz_KS, work3d_l, diag)
  endif
  if (VBF%id_Bdif_idz_KS>0) call post_data(VBF%id_Bdif_idz_KS, work2d, diag)
  if (VBF%id_Bdif_idV_KS>0) call post_data(VBF%id_Bdif_idV_KS, work, diag)

  ! Compute bkgnd fluxes
  if (VBF%id_Bdif_dz_bkgnd>0.or.VBF%id_Bdif_idz_bkgnd>0.or.VBF%id_Bdif_idV_bkgnd>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_bkgnd, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_bkgnd, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_idz_bkgnd>0) then
      work2d(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d(i,j) = work2d(i,j) + (VBF%Bflx_salt_dz(i,j,k) + VBF%Bflx_temp_dz(i,j,k))
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_bkgnd>0) then
      work = 0.0
      do k = 1,nz
        work = work + &
               (global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3) + &
                global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3))
      enddo
    endif
  elseif (VBF%id_Bdif_bkgnd>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_bkgnd, VBF%Bflx_salt)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_bkgnd, VBF%Bflx_temp)
  endif
  ! Post bkgnd fluxes
  if (VBF%id_Bdif_bkgnd>0) then
    work3d_i(:,:,:) = 0.0
    do k = 1,nz+1 ; do j = jsc,jec ; do i = isc,iec
      work3d_i(i,j,k) = VBF%Bflx_temp(i,j,k) + VBF%Bflx_salt(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_bkgnd, work3d_i, diag)
  endif
  if (VBF%id_Bdif_dz_bkgnd>0) then
    work3d_l(:,:,:) = 0.0
    do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
      work3d_l(i,j,k) = VBF%Bflx_temp_dz(i,j,k) + VBF%Bflx_salt_dz(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_dz_bkgnd, work3d_l, diag)
  endif
  if (VBF%id_Bdif_idz_bkgnd>0) call post_data(VBF%id_Bdif_idz_bkgnd, work2d, diag)
  if (VBF%id_Bdif_idV_bkgnd>0) call post_data(VBF%id_Bdif_idV_bkgnd, work, diag)

  ! Compute double diffusion fluxes
  if (VBF%id_Bdif_dz_ddiff_temp>0.or.VBF%id_Bdif_idz_ddiff_temp>0.or.VBF%id_Bdif_idV_ddiff_temp>0) then
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_ddiff_T, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_idz_ddiff_temp>0) then
      work2d_temp(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d_temp(i,j) = work2d_temp(i,j) + VBF%Bflx_temp_dz(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_ddiff_temp>0) then
      work_temp = 0.0
      do k = 1,nz
        work_temp = work_temp + global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, &
                    tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
      enddo
    endif
  elseif (VBF%id_Bdif_ddiff_temp>0) then
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_ddiff_T, VBF%Bflx_temp)
  endif
  if (VBF%id_Bdif_dz_ddiff_salt>0.or.VBF%id_Bdif_idz_ddiff_salt>0.or.VBF%id_Bdif_idV_ddiff_salt>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_ddiff_S, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    if (VBF%id_Bdif_idz_ddiff_salt>0) then
      work2d_salt(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d_salt(i,j) = work2d_salt(i,j) + VBF%Bflx_salt_dz(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_ddiff_salt>0) then
      work_salt = 0.0
      do k = 1,nz
        work_salt = work_salt + global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, &
                    tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
      enddo
    endif
  elseif (VBF%id_Bdif_ddiff_salt>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_ddiff_S, VBF%Bflx_salt)
  endif
  ! Post double diffusion fluxes
  if (VBF%id_Bdif_ddiff_temp>0) call post_data(VBF%id_Bdif_ddiff_temp, VBF%Bflx_temp, diag)
  if (VBF%id_Bdif_dz_ddiff_temp>0) call post_data(VBF%id_Bdif_dz_ddiff_temp, VBF%Bflx_temp_dz, diag)
  if (VBF%id_Bdif_idz_ddiff_temp>0) call post_data(VBF%id_Bdif_idz_ddiff_temp, work2d_temp, diag)
  if (VBF%id_Bdif_idV_ddiff_temp>0) call post_data(VBF%id_Bdif_idV_ddiff_temp, work_temp, diag)
  if (VBF%id_Bdif_ddiff_salt>0) call post_data(VBF%id_Bdif_ddiff_salt, VBF%Bflx_salt, diag)
  if (VBF%id_Bdif_dz_ddiff_salt>0) call post_data(VBF%id_Bdif_dz_ddiff_salt, VBF%Bflx_salt_dz, diag)
  if (VBF%id_Bdif_idz_ddiff_salt>0) call post_data(VBF%id_Bdif_idz_ddiff_salt, work2d_salt, diag)
  if (VBF%id_Bdif_idV_ddiff_salt>0) call post_data(VBF%id_Bdif_idV_ddiff_salt, work_salt, diag)

  ! Compute Kd_leak fluxes
  if (VBF%id_Bdif_dz_leak>0.or.VBF%id_Bdif_idz_leak>0.or.VBF%id_Bdif_idV_leak>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_leak, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_leak, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_idz_leak>0) then
      work2d(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d(i,j) = work2d(i,j) + (VBF%Bflx_salt_dz(i,j,k) + VBF%Bflx_temp_dz(i,j,k))
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_leak>0) then
      work = 0.0
      do k = 1,nz
        work = work + &
               (global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3) + &
                global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3))
      enddo
    endif
  elseif (VBF%id_Bdif_leak>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_leak, VBF%Bflx_salt)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_leak, VBF%Bflx_temp)
  endif
  ! Post Kd_leak fluxes
  if (VBF%id_Bdif_leak>0) then
    work3d_i(:,:,:) = 0.0
    do k = 1,nz+1 ; do j = jsc,jec ; do i = isc,iec
      work3d_i(i,j,k) = VBF%Bflx_temp(i,j,k) + VBF%Bflx_salt(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_leak, work3d_i, diag)
  endif
  if (VBF%id_Bdif_dz_leak>0) then
    work3d_l(:,:,:) = 0.0
    do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
      work3d_l(i,j,k) = VBF%Bflx_temp_dz(i,j,k) + VBF%Bflx_salt_dz(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_dz_leak, work3d_l, diag)
  endif
  if (VBF%id_Bdif_idz_leak>0) call post_data(VBF%id_Bdif_idz_leak, work2d, diag)
  if (VBF%id_Bdif_idV_leak>0) call post_data(VBF%id_Bdif_idV_leak, work, diag)

  ! Compute Kd_quad fluxes
  if (VBF%id_Bdif_dz_quad>0.or.VBF%id_Bdif_idz_quad>0.or.VBF%id_Bdif_idV_quad>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_quad, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_quad, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_idz_quad>0) then
      work2d(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d(i,j) = work2d(i,j) + (VBF%Bflx_salt_dz(i,j,k) + VBF%Bflx_temp_dz(i,j,k))
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_quad>0) then
      work = 0.0
      do k = 1,nz
        work = work + &
               (global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3) + &
                global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3))
      enddo
    endif
  elseif (VBF%id_Bdif_quad>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_quad, VBF%Bflx_salt)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_quad, VBF%Bflx_temp)
  endif
  ! Post Kd_quad fluxes
  if (VBF%id_Bdif_quad>0) then
    work3d_i(:,:,:) = 0.0
    do k = 1,nz+1 ; do j = jsc,jec ; do i = isc,iec
      work3d_i(i,j,k) = VBF%Bflx_temp(i,j,k) + VBF%Bflx_salt(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_quad, work3d_i, diag)
  endif
  if (VBF%id_Bdif_dz_quad>0) then
    work3d_l(:,:,:) = 0.0
    do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
      work3d_l(i,j,k) = VBF%Bflx_temp_dz(i,j,k) + VBF%Bflx_salt_dz(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_dz_quad, work3d_l, diag)
  endif
  if (VBF%id_Bdif_idz_quad>0) call post_data(VBF%id_Bdif_idz_quad, work2d, diag)
  if (VBF%id_Bdif_idV_quad>0) call post_data(VBF%id_Bdif_idV_quad, work, diag)

  ! Compute Kd_itidal fluxes
  if (VBF%id_Bdif_dz_itidal>0.or.VBF%id_Bdif_idz_itidal>0.or.VBF%id_Bdif_idV_itidal>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_itidal, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_itidal, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_idz_itidal>0) then
      work2d(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d(i,j) = work2d(i,j) + (VBF%Bflx_salt_dz(i,j,k) + VBF%Bflx_temp_dz(i,j,k))
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_itidal>0) then
      work = 0.0
      do k = 1,nz
        work = work + &
               (global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3) + &
                global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3))
      enddo
    endif
  elseif (VBF%id_Bdif_itidal>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_itidal, VBF%Bflx_salt)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_itidal, VBF%Bflx_temp)
  endif
  ! Post Kd_itidal fluxes
  if (VBF%id_Bdif_itidal>0) then
    work3d_i(:,:,:) = 0.0
    do k = 1,nz+1 ; do j = jsc,jec ; do i = isc,iec
      work3d_i(i,j,k) = VBF%Bflx_temp(i,j,k) + VBF%Bflx_salt(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_itidal, work3d_i, diag)
  endif
  if (VBF%id_Bdif_dz_itidal>0) then
    work3d_l(:,:,:) = 0.0
    do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
      work3d_l(i,j,k) = VBF%Bflx_temp_dz(i,j,k)+VBF%Bflx_salt_dz(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_dz_itidal, work3d_l, diag)
  endif
  if (VBF%id_Bdif_idz_itidal>0) call post_data(VBF%id_Bdif_idz_itidal, work2d, diag)
  if (VBF%id_Bdif_idV_itidal>0) call post_data(VBF%id_Bdif_idV_itidal, work, diag)

  ! Compute Kd_Froude fluxes
  if (VBF%id_Bdif_dz_Froude>0.or.VBF%id_Bdif_idz_Froude>0.or.VBF%id_Bdif_idV_Froude>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_Froude, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_Froude, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_idz_Froude>0) then
      work2d(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d(i,j) = work2d(i,j) + (VBF%Bflx_salt_dz(i,j,k) + VBF%Bflx_temp_dz(i,j,k))
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_Froude>0) then
      work = 0.0
      do k = 1,nz
        work = work + &
               (global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3) + &
                global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3))
      enddo
    endif
  elseif (VBF%id_Bdif_Froude>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_Froude, VBF%Bflx_salt)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_Froude, VBF%Bflx_temp)
  endif
  ! Post Kd_Froude fluxes
  if (VBF%id_Bdif_Froude>0) then
    work3d_i(:,:,:) = 0.0
    do k = 1,nz+1 ; do j = jsc,jec ; do i = isc,iec
      work3d_i(i,j,k) = VBF%Bflx_temp(i,j,k) + VBF%Bflx_salt(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_Froude, work3d_i, diag)
  endif
  if (VBF%id_Bdif_dz_Froude>0) then
    work3d_l(:,:,:) = 0.0
    do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
      work3d_l(i,j,k) = VBF%Bflx_temp_dz(i,j,k) + VBF%Bflx_salt_dz(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_dz_Froude, work3d_l, diag)
  endif
  if (VBF%id_Bdif_idz_Froude>0) call post_data(VBF%id_Bdif_idz_Froude, work2d, diag)
  if (VBF%id_Bdif_idV_Froude>0) call post_data(VBF%id_Bdif_idV_Froude, work, diag)

  ! Compute Kd_slope fluxes
  if (VBF%id_Bdif_dz_slope>0.or.VBF%id_Bdif_idz_slope>0.or.VBF%id_Bdif_idV_slope>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_slope, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_slope, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_idz_slope>0) then
      work2d(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d(i,j) = work2d(i,j) + (VBF%Bflx_salt_dz(i,j,k) + VBF%Bflx_temp_dz(i,j,k))
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_slope>0) then
      work = 0.0
      do k = 1,nz
        work = work + &
               (global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3) + &
                global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3))
      enddo
    endif
  elseif (VBF%id_Bdif_slope>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_slope, VBF%Bflx_salt)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_slope, VBF%Bflx_temp)
  endif
  ! Post Kd_slope fluxes
  if (VBF%id_Bdif_slope>0) then
    work3d_i(:,:,:) = 0.0
    do k = 1,nz+1 ; do j = jsc,jec ; do i = isc,iec
      work3d_i(i,j,k) = VBF%Bflx_temp(i,j,k) + VBF%Bflx_salt(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_slope, work3d_i, diag)
  endif
  if (VBF%id_Bdif_dz_slope>0) then
    work3d_l(:,:,:) = 0.0
    do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
      work3d_l(i,j,k) = VBF%Bflx_temp_dz(i,j,k) + VBF%Bflx_salt_dz(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_dz_slope, work3d_l, diag)
  endif
  if (VBF%id_Bdif_idz_slope>0) call post_data(VBF%id_Bdif_idz_slope, work2d, diag)
  if (VBF%id_Bdif_idV_slope>0) call post_data(VBF%id_Bdif_idV_slope, work, diag)

  ! Compute Kd_lowmode fluxes
  if (VBF%id_Bdif_dz_lowmode>0.or.VBF%id_Bdif_idz_lowmode>0.or.VBF%id_Bdif_idV_lowmode>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_lowmode, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_lowmode, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_idz_lowmode>0) then
      work2d(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d(i,j) = work2d(i,j) + (VBF%Bflx_salt_dz(i,j,k) + VBF%Bflx_temp_dz(i,j,k))
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_lowmode>0) then
      work = 0.0
      do k = 1,nz
        work = work + &
               (global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3) + &
                global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3))
      enddo
    endif
  elseif (VBF%id_Bdif_lowmode>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_lowmode, VBF%Bflx_salt)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_lowmode, VBF%Bflx_temp)
  endif
  ! Post Kd_lowmode fluxes
  if (VBF%id_Bdif_lowmode>0) then
    work3d_i(:,:,:) = 0.0
    do k = 1,nz+1 ; do j = jsc,jec ; do i = isc,iec
      work3d_i(i,j,k) = VBF%Bflx_temp(i,j,k) + VBF%Bflx_salt(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_lowmode, work3d_i, diag)
  endif
  if (VBF%id_Bdif_dz_lowmode>0) then
    work3d_l(:,:,:) = 0.0
    do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
      work3d_l(i,j,k) = VBF%Bflx_temp_dz(i,j,k) + VBF%Bflx_salt_dz(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_dz_lowmode, work3d_l, diag)
  endif
  if (VBF%id_Bdif_idz_lowmode>0) call post_data(VBF%id_Bdif_idz_lowmode, work2d, diag)
  if (VBF%id_Bdif_idV_lowmode>0) call post_data(VBF%id_Bdif_idV_lowmode, work, diag)

  ! Compute Kd_Niku fluxes
  if (VBF%id_Bdif_dz_Niku>0 .or. VBF%id_Bdif_idz_Niku>0 .or. VBF%id_Bdif_idV_Niku>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_Niku, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_Niku, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_idz_Niku>0) then
      work2d(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d(i,j) = work2d(i,j) + (VBF%Bflx_salt_dz(i,j,k) + VBF%Bflx_temp_dz(i,j,k))
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_Niku>0) then
      work = 0.0
      do k = 1,nz
        work = work + &
               (global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3) + &
                global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3))
      enddo
    endif
  elseif (VBF%id_Bdif_Niku>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_Niku, VBF%Bflx_salt)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_Niku, VBF%Bflx_temp)
  endif
  ! Post Kd_Niku fluxes
  if (VBF%id_Bdif_Niku>0) then
    work3d_i(:,:,:) = 0.0
    do k = 1,nz+1 ; do j = jsc,jec ; do i = isc,iec
      work3d_i(i,j,k) = VBF%Bflx_temp(i,j,k) + VBF%Bflx_salt(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_lowmode, work3d_i, diag)
  endif
  if (VBF%id_Bdif_dz_Niku>0) then
    work3d_l(:,:,:) = 0.0
    do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
      work3d_l(i,j,k) = VBF%Bflx_temp_dz(i,j,k) + VBF%Bflx_salt_dz(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_dz_Niku, work3d_l, diag)
  endif
  if (VBF%id_Bdif_idz_Niku>0) call post_data(VBF%id_Bdif_idz_Niku, work2d, diag)
  if (VBF%id_Bdif_idV_Niku>0) call post_data(VBF%id_Bdif_idV_Niku, work, diag)

  ! Compute Kd_itides fluxes
  if (VBF%id_Bdif_dz_itides>0 .or. VBF%id_Bdif_idz_itides>0 .or. VBF%id_Bdif_idV_itides>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_itides, VBF%Bflx_salt, dz=dz, Bdif_flx_dz=VBF%Bflx_salt_dz)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_itides, VBF%Bflx_temp, dz=dz, Bdif_flx_dz=VBF%Bflx_temp_dz)
    if (VBF%id_Bdif_idz_itides>0) then
      work2d(:,:) = 0.0
      do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
        work2d(i,j) = work2d(i,j) + (VBF%Bflx_salt_dz(i,j,k) + VBF%Bflx_temp_dz(i,j,k))
      enddo ; enddo ; enddo
    endif
    if (VBF%id_Bdif_idV_itides>0) then
      work = 0.0
      do k = 1,nz
        work = work + &
               (global_area_integral(VBF%Bflx_temp_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3) + &
                global_area_integral(VBF%Bflx_salt_dz(:,:,k), G, tmp_scale=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3))
      enddo
    endif
  elseif (VBF%id_Bdif_itides>0) then
    call diagnoseKdWork(G, GV, N2_salt, VBF%Kd_itides, VBF%Bflx_salt)
    call diagnoseKdWork(G, GV, N2_temp, VBF%Kd_itides, VBF%Bflx_temp)
  endif
  ! Post Kd_itides fluxes
  if (VBF%id_Bdif_itides>0) then
    work3d_i(:,:,:) = 0.0
    do k = 1,nz+1 ; do j = jsc,jec ; do i = isc,iec
      work3d_i(i,j,k) = VBF%Bflx_temp(i,j,k) + VBF%Bflx_salt(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_itides, work3d_i, diag)
  endif
  if (VBF%id_Bdif_dz_itides>0) then
    work3d_l(:,:,:) = 0.0
    do k = 1,nz ; do j = jsc,jec ; do i = isc,iec
      work3d_l(i,j,k) = VBF%Bflx_temp_dz(i,j,k) + VBF%Bflx_salt_dz(i,j,k)
    enddo ; enddo ; enddo
    call post_data(VBF%id_Bdif_dz_itides, work3d_l, diag)
  endif
  if (VBF%id_Bdif_idz_itides>0) call post_data(VBF%id_Bdif_idz_itides, work2d, diag)
  if (VBF%id_Bdif_idV_itides>0) call post_data(VBF%id_Bdif_idV_itides, work, diag)

end subroutine KdWork_Diagnostics

!> Diagnose the implied "work", or buoyancy forcing & its integral, due to a given diffusivity and column state.
subroutine diagnoseKdWork(G, GV, N2, Kd, Bdif_flx, dz, Bdif_flx_dz)
  type(ocean_grid_type),   intent(in)  :: G    !< Grid type
  type(verticalGrid_type), intent(in)  :: GV   !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(in)  :: N2   !< Buoyancy frequency [T-2 ~> s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(in)  :: Kd   !< Diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(out) :: Bdif_flx !< Buoyancy flux [H Z T-3 ~> m2 s-3 or kg m-1 s-3 = W m-3]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 intent(in), optional :: dz    !< Grid spacing [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                intent(out), optional :: Bdif_flx_dz !< Buoyancy flux over layer [H Z2 T-3 ~> m3 s-3 or kg s-3 = W m-2]

  integer :: i, j, k

  Bdif_flx(:,:,1) = 0.0
  Bdif_flx(:,:,GV%ke+1) = 0.0
  !$OMP parallel do default(shared)
  do K=2,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    Bdif_flx(i,j,K) = - N2(i,j,K) * Kd(i,j,K)
  enddo ; enddo; enddo

  if (present(Bdif_flx_dz) .and. present(dz)) then
    !$OMP parallel do default(shared)
    do K=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      Bdif_flx_dz(i,j,k) = 0.5*(Bdif_flx(i,j,K)+Bdif_flx(i,j,K+1))*dz(i,j,k)
    enddo ; enddo; enddo
  endif

end subroutine diagnoseKdWork

!> Allocates arrays only when needed
subroutine Allocate_VBF_CS(G, GV, VBF)
  type(ocean_grid_type),   intent(in) :: G   !< ocean grid structure
  type(verticalGrid_type), intent(in) :: GV  !< ocean vertical grid structure
  type (vbf_CS),        intent(inout) :: VBF !< Vertical buoyancy flux structure

  integer :: isd, ied, jsd, jed, nz

  isd  = G%isd  ; ied = G%ied  ; jsd = G%jsd  ; jed = G%jed ; nz = GV%ke

  if (VBF%do_bflx_salt) &
    allocate(VBF%Bflx_salt(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%do_bflx_salt_dz) &
     allocate(VBF%Bflx_salt_dz(isd:ied,jsd:jed,nz), source=0.0)
  if (VBF%do_bflx_temp) &
     allocate(VBF%Bflx_temp(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%do_bflx_temp_dz) &
     allocate(VBF%Bflx_temp_dz(isd:ied,jsd:jed,nz), source=0.0)

  if (VBF%id_Bdif_salt_dz>0 .or. VBF%id_Bdif_dz>0 .or. VBF%id_Bdif_salt>0 .or. VBF%id_Bdif>0 .or. &
      VBF%id_Bdif_idz>0 .or. VBF%id_Bdif_salt_idz>0 .or. VBF%id_Bdif_idV>0 .or. VBF%id_Bdif_salt_idV>0) &
    allocate(VBF%Kd_salt(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_temp_dz>0 .or. VBF%id_Bdif_dz>0 .or. VBF%id_Bdif_temp>0 .or. VBF%id_Bdif>0 .or. &
      VBF%id_Bdif_idz>0 .or. VBF%id_Bdif_temp_idz>0 .or. VBF%id_Bdif_idV>0 .or. VBF%id_Bdif_temp_idV>0) &
    allocate(VBF%Kd_temp(isd:ied,jsd:jed,nz+1), source=0.0)

  if (VBF%id_Bdif_BBL>0 .or. VBF%id_Bdif_dz_BBL>0 .or. VBF%id_Bdif_idz_BBL>0 .or. VBF%id_Bdif_idV_BBL>0) &
    allocate(VBF%Kd_BBL(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_ePBL>0 .or. VBF%id_Bdif_dz_ePBL>0 .or. VBF%id_Bdif_idz_ePBL>0 .or. VBF%id_Bdif_idV_ePBL>0) &
    allocate(VBF%Kd_ePBL(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_KS>0 .or. VBF%id_Bdif_dz_KS>0 .or. VBF%id_Bdif_idz_KS>0 .or. VBF%id_Bdif_idV_KS>0) &
    allocate(VBF%Kd_KS(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_bkgnd>0 .or. VBF%id_Bdif_dz_bkgnd>0 .or. VBF%id_Bdif_idz_bkgnd>0 .or. VBF%id_Bdif_idV_bkgnd>0) &
    allocate(VBF%Kd_bkgnd(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_ddiff_temp>0 .or. VBF%id_Bdif_dz_ddiff_temp>0 .or. VBF%id_Bdif_idz_ddiff_temp>0 &
      .or. VBF%id_Bdif_idV_ddiff_temp>0) allocate(VBF%Kd_ddiff_T(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_ddiff_salt>0 .or. VBF%id_Bdif_dz_ddiff_salt>0 .or. VBF%id_Bdif_idV_ddiff_salt>0 &
      .or. VBF%id_Bdif_idV_ddiff_salt>0) allocate(VBF%Kd_ddiff_S(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_leak>0 .or. VBF%id_Bdif_dz_leak>0 .or. VBF%id_Bdif_idz_leak>0 .or. VBF%id_Bdif_idV_leak>0) &
    allocate(VBF%Kd_leak(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_quad>0 .or. VBF%id_Bdif_dz_quad>0 .or. VBF%id_Bdif_idz_quad>0 .or. VBF%id_Bdif_idV_quad>0) &
    allocate(VBF%Kd_quad(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_itidal>0 .or. VBF%id_Bdif_dz_itidal>0 .or. VBF%id_Bdif_idz_itidal>0 .or. VBF%id_Bdif_idV_itidal>0) &
    allocate(VBF%Kd_itidal(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_Froude>0 .or. VBF%id_Bdif_dz_Froude>0 .or. VBF%id_Bdif_idz_Froude>0 .or. VBF%id_Bdif_idV_Froude>0) &
    allocate(VBF%Kd_Froude(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_slope>0 .or. VBF%id_Bdif_dz_slope>0 .or. VBF%id_Bdif_idz_slope>0 .or. VBF%id_Bdif_idV_slope>0) &
    allocate(VBF%Kd_slope(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_lowmode>0 .or. VBF%id_Bdif_dz_lowmode>0 .or. VBF%id_Bdif_idz_lowmode>0 .or. &
      VBF%id_Bdif_idV_lowmode>0) allocate(VBF%Kd_lowmode(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_Niku>0 .or. VBF%id_Bdif_dz_Niku>0 .or. VBF%id_Bdif_idz_Niku>0 .or. VBF%id_Bdif_idV_Niku>0) &
    allocate(VBF%Kd_Niku(isd:ied,jsd:jed,nz+1), source=0.0)
  if (VBF%id_Bdif_itides>0 .or. VBF%id_Bdif_dz_itides>0 .or. VBF%id_Bdif_idz_itides>0 .or. VBF%id_Bdif_idV_itides>0) &
    allocate(VBF%Kd_itides(isd:ied,jsd:jed,nz+1), source=0.0)

end subroutine Allocate_VBF_CS

!> Deallocate any arrays that were allocated
subroutine Deallocate_VBF_CS(VBF)
  type (vbf_CS), intent(inout) :: VBF !< Vertical buoyancy flux structure

  if (associated(VBF%Bflx_salt)) &
    deallocate(VBF%Bflx_salt)
  if (associated(VBF%Bflx_temp)) &
    deallocate(VBF%Bflx_temp)
  if (associated(VBF%Bflx_salt_dz)) &
    deallocate(VBF%Bflx_salt_dz)
  if (associated(VBF%Bflx_temp_dz)) &
    deallocate(VBF%Bflx_temp_dz)
  if (associated(VBF%Kd_salt)) &
    deallocate(VBF%Kd_salt)
  if (associated(VBF%Kd_temp)) &
    deallocate(VBF%Kd_temp)
  if (associated(VBF%Kd_BBL)) &
    deallocate(VBF%Kd_BBL)
  if (associated(VBF%Kd_ePBL)) &
    deallocate(VBF%Kd_ePBL)
  if (associated(VBF%Kd_KS)) &
    deallocate(VBF%Kd_KS)
  if (associated(VBF%Kd_bkgnd)) &
    deallocate(VBF%Kd_bkgnd)
  if (associated(VBF%Kd_ddiff_T)) &
    deallocate(VBF%Kd_ddiff_T)
  if (associated(VBF%Kd_ddiff_S)) &
    deallocate(VBF%Kd_ddiff_S)
  if (associated(VBF%Kd_leak)) &
    deallocate(VBF%Kd_leak)
  if (associated(VBF%Kd_quad)) &
    deallocate(VBF%Kd_quad)
  if (associated(VBF%Kd_itidal)) &
    deallocate(VBF%Kd_itidal)
  if (associated(VBF%Kd_Froude)) &
    deallocate(VBF%Kd_Froude)
  if (associated(VBF%Kd_slope)) &
    deallocate(VBF%Kd_slope)
  if (associated(VBF%Kd_lowmode)) &
    deallocate(VBF%Kd_lowmode)
  if (associated(VBF%Kd_Niku)) &
    deallocate(VBF%Kd_Niku)
  if (associated(VBF%Kd_itides)) &
    deallocate(VBF%Kd_itides)

end subroutine Deallocate_VBF_CS

!> Handles all KdWork diagnostics and flags which calculations should be done.
subroutine KdWork_init(Time, G,GV,US,diag,VBF,Use_KdWork_diag)
  type(time_type), target                :: Time             !< model time
  type(ocean_grid_type),   intent(in)    :: G        !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV       !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US       !< A dimensional unit scaling type
  type(diag_ctrl), target, intent(inout) :: diag     !< regulates diagnostic output
  type (vbf_CS), pointer,  intent(inout) :: VBF      !< Vertical buoyancy flux structure
  logical,                 intent(out)   :: Use_KdWork_diag !< Flag if any output was turned on

  allocate(VBF)

  VBF%do_bflx_salt = .false.
  VBF%do_bflx_salt_dz = .false.
  VBF%do_bflx_temp = .false.
  VBF%do_bflx_temp_dz = .false.

  VBF%id_Bdif = register_diag_field('ocean_model',"Bflx_dia_diff", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz = register_diag_field('ocean_model',"Bflx_dia_diff_dz", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz = register_diag_field('ocean_model',"Bflx_dia_diff_idz", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV = register_scalar_field('ocean_model',"Bflx_dia_diff_idV", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_salt = register_diag_field('ocean_model',"Bflx_salt_dia_diff", diag%axesTi, &
        Time, "Salinity contribution to diffusive diapycnal buoyancy flux across interfaces", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_salt_dz = register_diag_field('ocean_model',"Bflx_salt_dia_diff_dz", diag%axesTl, &
        Time, "Salinity contribution to layer integral of diffusive diapycnal buoyancy flux.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_salt_idz = register_diag_field('ocean_model',"Bflx_salt_dia_diff_idz", diag%axesT1, &
        Time, "Salinity contribution to layer integrated diffusive diapycnal buoyancy flux.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_salt_idV = register_scalar_field('ocean_model',"Bflx_salt_dia_diff_idV", Time, diag, &
        "Salinity contribution to global integrated diffusive diapycnal buoyancy flux.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_temp = register_diag_field('ocean_model',"Bflx_temp_dia_diff", diag%axesTi, &
        Time, "Temperature contribution to diffusive diapycnal buoyancy flux across interfaces", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_temp_dz = register_diag_field('ocean_model',"Bflx_temp_dia_diff_dz", diag%axesTl, &
        Time, "Temperature contribution to layer integral of diffusive diapycnal buoyancy flux.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_temp_idz = register_diag_field('ocean_model',"Bflx_temp_dia_diff_idz", diag%axesT1, &
        Time, "Temperature contribution to layer integrated diffusive diapycnal buoyancy flux.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_temp_idV = register_scalar_field('ocean_model',"Bflx_temp_dia_diff_idV", Time, diag, &
        "Temperature contribution to global integrated diffusive diapycnal buoyancy flux.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_BBL = register_diag_field('ocean_model',"Bflx_dia_diff_BBL", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to the BBL parameterization.", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_BBL = register_diag_field('ocean_model',"Bflx_dia_diff_dz_BBL", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to the BBL parameterization.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_BBL = register_diag_field('ocean_model',"Bflx_dia_diff_idz_BBL", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to the BBL parameterization.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_BBL = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_BBL", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to BBL.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_ePBL = register_diag_field('ocean_model',"Bflx_dia_diff_ePBL", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to ePBL", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_ePBL = register_diag_field('ocean_model',"Bflx_dia_diff_dz_ePBL", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to ePBL.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_ePBL = register_diag_field('ocean_model',"Bflx_dia_diff_idz_ePBL", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to ePBL.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_ePBL = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_ePBL", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to ePBL.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_KS = register_diag_field('ocean_model',"Bflx_dia_diff_KS", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to Kappa Shear", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_KS = register_diag_field('ocean_model',"Bflx_dia_diff_dz_KS", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to Kappa Shear.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_KS = register_diag_field('ocean_model',"Bflx_dia_diff_idz_KS", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to Kappa Shear.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_KS = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_KS", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to Kappa Shear.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_bkgnd = register_diag_field('ocean_model',"Bflx_dia_diff_bkgnd", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to bkgnd mixing", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_bkgnd = register_diag_field('ocean_model',"Bflx_dia_diff_dz_bkgnd", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_bkgnd = register_diag_field('ocean_model',"Bflx_dia_diff_idz_bkgnd", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_bkgnd = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_bkgnd", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to Kd_bkgnd.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_ddiff_temp = register_diag_field('ocean_model',"Bflx_dia_diff_ddiff_heat", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to double diffusion of heat", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_ddiff_temp = register_diag_field('ocean_model',"Bflx_dia_diff_dz_ddiff_heat", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to double diffusion of heat.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_ddiff_temp = register_diag_field('ocean_model',"Bflx_dia_diff_idz_ddiff_heat", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to double diffusion of heat.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_ddiff_temp = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_ddiff_heat", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to double diffusion of heat.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_ddiff_salt = register_diag_field('ocean_model',"Bflx_dia_diff_ddiff_salt", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to double diffusion of salt", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_ddiff_salt = register_diag_field('ocean_model',"Bflx_dia_diff_dz_ddiff_salt", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to double diffusion of salt.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_ddiff_salt = register_diag_field('ocean_model',"Bflx_dia_diff_idz_ddiff_salt", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to double diffusion of salt.", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_ddiff_salt = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_ddiff_salt", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to double diffusion of salt.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_leak = register_diag_field('ocean_model',"Bflx_dia_diff_leak", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to Kd_leak mixing", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_leak = register_diag_field('ocean_model',"Bflx_dia_diff_dz_leak", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_leak = register_diag_field('ocean_model',"Bflx_dia_diff_idz_leak", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_leak = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_leak", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to Kd_leak.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_quad = register_diag_field('ocean_model',"Bflx_dia_diff_quad", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to Kd_quad mixing", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_quad = register_diag_field('ocean_model',"Bflx_dia_diff_dz_quad", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_quad = register_diag_field('ocean_model',"Bflx_dia_diff_idz_quad", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_quad = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_quad", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to Kd_quad.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_itidal = register_diag_field('ocean_model',"Bflx_dia_diff_itidal", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to Kd_itidal mixing", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_itidal = register_diag_field('ocean_model',"Bflx_dia_diff_dz_itidal", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_itidal = register_diag_field('ocean_model',"Bflx_dia_diff_idz_itidal", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_itidal = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_itidal", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to Kd_itidal.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_Froude = register_diag_field('ocean_model',"Bflx_dia_diff_Froude", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to Kd_Froude mixing", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_Froude = register_diag_field('ocean_model',"Bflx_dia_diff_dz_Froude", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_Froude = register_diag_field('ocean_model',"Bflx_dia_diff_idz_Froude", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_Froude = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_Froude", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to Kd_Froude.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_slope = register_diag_field('ocean_model',"Bflx_dia_diff_slope", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to Kd_slope mixing", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_slope = register_diag_field('ocean_model',"Bflx_dia_diff_dz_slope", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_slope = register_diag_field('ocean_model',"Bflx_dia_diff_idz_slope", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_slope = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_slope", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to Kd_slope.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_lowmode = register_diag_field('ocean_model',"Bflx_dia_diff_lowmode", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to Kd_lowmode mixing", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_lowmode = register_diag_field('ocean_model',"Bflx_dia_diff_dz_lowmode", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_lowmode = register_diag_field('ocean_model',"Bflx_dia_diff_idz_lowmode", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_lowmode = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_lowmode", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to Kd_lowmode.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_Niku = register_diag_field('ocean_model',"Bflx_dia_diff_Niku", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to Kd_Niku mixing", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_Niku = register_diag_field('ocean_model',"Bflx_dia_diff_dz_Niku", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_Niku = register_diag_field('ocean_model',"Bflx_dia_diff_idz_Niku", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_Niku = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_Niku", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to Kd_Niku.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  VBF%id_Bdif_itides = register_diag_field('ocean_model',"Bflx_dia_diff_itides", diag%axesTi, &
        Time, "Diffusive diapycnal buoyancy flux across interfaces due to Kd_itides mixing", &
        "W m-3", conversion=GV%H_to_kg_m2*US%Z_to_m*US%s_to_T**3)
  VBF%id_Bdif_dz_itides = register_diag_field('ocean_model',"Bflx_dia_diff_dz_itides", diag%axesTl, &
        Time, "Layerwise integral of diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idz_itides = register_diag_field('ocean_model',"Bflx_dia_diff_idz_itides", diag%axesT1, &
        Time, "Layer integrated diffusive diapycnal buoyancy flux due to bkgnd mixing", &
        "W m-2", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3)
  VBF%id_Bdif_idV_itides = register_scalar_field('ocean_model',"Bflx_dia_diff_idV_itides", Time, diag, &
        "Global integrated diffusive diapycnal buoyancy flux due to Kd_itides.", &
        units="W", conversion=GV%H_to_kg_m2*US%Z_to_m**2*US%s_to_T**3*US%L_to_m**2)

  if (VBF%id_Bdif_dz>0 .or. VBF%id_Bdif_salt_dz>0 .or. VBF%id_Bdif_dz_BBL>0 .or. &
      VBF%id_Bdif_dz_ePBL>0 .or. VBF%id_Bdif_dz_KS>0 .or. VBF%id_Bdif_dz_bkgnd>0 .or. &
      VBF%id_Bdif_dz_ddiff_salt>0 .or. VBF%id_Bdif_dz_leak>0 .or. VBF%id_Bdif_dz_quad>0 .or. &
      VBF%id_Bdif_dz_itidal>0 .or. VBF%id_Bdif_dz_Froude>0 .or. VBF%id_Bdif_dz_slope>0 .or. &
      VBF%id_Bdif_dz_lowmode>0 .or. VBF%id_Bdif_dz_Niku>0 .or. VBF%id_Bdif_dz_itides>0 .or. &
      VBF%id_Bdif_idV>0 .or. VBF%id_Bdif_salt_idV>0 .or. VBF%id_Bdif_idV_BBL>0 .or. &
      VBF%id_Bdif_idV_ePBL>0 .or. VBF%id_Bdif_idV_KS>0 .or. VBF%id_Bdif_idV_bkgnd>0 .or. &
      VBF%id_Bdif_idV_ddiff_salt>0 .or. VBF%id_Bdif_idV_leak>0 .or. VBF%id_Bdif_idV_quad>0 .or. &
      VBF%id_Bdif_idV_itidal>0 .or. VBF%id_Bdif_idV_Froude>0 .or. VBF%id_Bdif_idV_slope>0 .or. &
      VBF%id_Bdif_idV_lowmode>0 .or. VBF%id_Bdif_idV_Niku>0 .or. VBF%id_Bdif_idV_itides>0 .or. &
      VBF%id_Bdif_idz>0 .or. VBF%id_Bdif_salt_idz>0 .or. VBF%id_Bdif_idz_BBL>0 .or. &
      VBF%id_Bdif_idz_ePBL>0 .or. VBF%id_Bdif_idz_KS>0 .or. VBF%id_Bdif_idz_bkgnd>0 .or. &
      VBF%id_Bdif_idz_ddiff_salt>0 .or. VBF%id_Bdif_idz_leak>0 .or. VBF%id_Bdif_idz_quad>0 .or. &
      VBF%id_Bdif_idz_itidal>0 .or. VBF%id_Bdif_idz_Froude>0 .or. VBF%id_Bdif_idz_slope>0 .or. &
      VBF%id_Bdif_idz_lowmode>0 .or. VBF%id_Bdif_idz_Niku>0 .or. VBF%id_Bdif_idz_itides>0 ) then
    VBF%do_bflx_salt_dz = .true.
  endif
  if (VBF%id_Bdif_dz>0 .or. VBF%id_Bdif_temp_dz>0 .or. VBF%id_Bdif_dz_BBL>0 .or. &
      VBF%id_Bdif_dz_ePBL>0 .or. VBF%id_Bdif_dz_KS>0 .or. VBF%id_Bdif_dz_bkgnd>0 .or. &
      VBF%id_Bdif_dz_ddiff_temp>0 .or. VBF%id_Bdif_dz_leak>0 .or. VBF%id_Bdif_dz_quad>0 .or. &
      VBF%id_Bdif_dz_itidal>0 .or. VBF%id_Bdif_dz_Froude>0 .or. VBF%id_Bdif_dz_slope>0 .or. &
      VBF%id_Bdif_dz_lowmode>0 .or. VBF%id_Bdif_dz_Niku>0 .or. VBF%id_Bdif_dz_itides>0 .or. &
      VBF%id_Bdif_idV>0 .or. VBF%id_Bdif_temp_idV>0 .or. VBF%id_Bdif_idV_BBL>0 .or. &
      VBF%id_Bdif_idV_ePBL>0 .or. VBF%id_Bdif_idV_KS>0 .or. VBF%id_Bdif_idV_bkgnd>0 .or. &
      VBF%id_Bdif_idV_ddiff_temp>0 .or. VBF%id_Bdif_idV_leak>0 .or. VBF%id_Bdif_idV_quad>0 .or. &
      VBF%id_Bdif_idV_itidal>0 .or. VBF%id_Bdif_idV_Froude>0 .or. VBF%id_Bdif_idV_slope>0 .or. &
      VBF%id_Bdif_idV_lowmode>0 .or. VBF%id_Bdif_idV_Niku>0 .or. VBF%id_Bdif_idV_itides>0 .or. &
      VBF%id_Bdif_idz>0 .or. VBF%id_Bdif_temp_idz>0 .or. VBF%id_Bdif_idz_BBL>0 .or. &
      VBF%id_Bdif_idz_ePBL>0 .or. VBF%id_Bdif_idz_KS>0 .or. VBF%id_Bdif_idz_bkgnd>0 .or. &
      VBF%id_Bdif_idz_ddiff_temp>0 .or. VBF%id_Bdif_idz_leak>0 .or. VBF%id_Bdif_idz_quad>0 .or. &
      VBF%id_Bdif_idz_itidal>0 .or. VBF%id_Bdif_idz_Froude>0 .or. VBF%id_Bdif_idz_slope>0 .or. &
      VBF%id_Bdif_idz_lowmode>0 .or. VBF%id_Bdif_idz_Niku>0 .or. VBF%id_Bdif_idz_itides>0 ) then
    VBF%do_bflx_temp_dz = .true.
  endif
  if (VBF%id_Bdif>0 .or. VBF%id_Bdif_salt>0 .or. VBF%id_Bdif_BBL>0 .or. &
      VBF%id_Bdif_ePBL>0 .or. VBF%id_Bdif_KS>0 .or. VBF%id_Bdif_bkgnd>0 .or. &
      VBF%id_Bdif_ddiff_salt>0 .or. VBF%id_Bdif_leak>0 .or. VBF%id_Bdif_quad>0 .or. &
      VBF%id_Bdif_itidal>0 .or. VBF%id_Bdif_Froude>0 .or. VBF%id_Bdif_slope>0 .or. &
      VBF%id_Bdif_lowmode>0 .or. VBF%id_Bdif_Niku>0 .or. VBF%id_Bdif_itides>0 .or. &
      VBF%do_bflx_salt_dz) then
    VBF%do_bflx_salt = .true.
  endif
  if (VBF%id_Bdif>0 .or. VBF%id_Bdif_temp>0 .or. VBF%id_Bdif_BBL>0 .or. &
      VBF%id_Bdif_ePBL>0 .or. VBF%id_Bdif_KS>0 .or. VBF%id_Bdif_bkgnd>0 .or. &
      VBF%id_Bdif_ddiff_temp>0 .or. VBF%id_Bdif_leak>0 .or. VBF%id_Bdif_quad>0 .or. &
      VBF%id_Bdif_itidal>0 .or. VBF%id_Bdif_Froude>0 .or. VBF%id_Bdif_slope>0 .or. &
      VBF%id_Bdif_lowmode>0 .or. VBF%id_Bdif_Niku>0 .or. VBF%id_Bdif_itides>0 .or. &
      VBF%do_bflx_temp_dz) then
    VBF%do_bflx_temp = .true.
  endif

  Use_KdWork_diag = (VBF%do_bflx_salt .or. VBF%do_bflx_temp .or. VBF%do_bflx_salt_dz .or. VBF%do_bflx_temp_dz)

end subroutine KdWork_init

!> Deallocates control structrue
subroutine KdWork_end(VBF)
  type (vbf_CS), pointer,  intent(inout) :: VBF      !< Vertical buoyancy flux structure

  if (associated(VBF)) deallocate(VBF)

end subroutine KdWork_end

!> \namespace mom_diagnose_kdwork
!!
!!    The subroutine diagnoseKdWork diagnoses the energetics associated with various vertical diffusivities
!!    inside MOM6 diabatic routines.
!!

end module MOM_diagnose_kdwork
