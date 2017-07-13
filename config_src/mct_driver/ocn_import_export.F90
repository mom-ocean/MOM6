module ocn_import_export

   implicit none
   public
   save

   ! accumulated sum of send buffer quantities for averaging before being sent
   !real (r8), dimension(:,:,:,:), allocatable ::  SBUFF_SUM 
   !real (r8) :: tlast_coupled 

   !TODO: update the types of following vars
   double precision, dimension(:,:,:,:), allocatable ::  SBUFF_SUM 
   double precision :: tlast_coupled 
contains

!***********************************************************************
!BOP
! !IROUTINE: ocn_import
! !INTERFACE:

  subroutine ocn_import(x2o, ldiag_cpl, errorCode)

! !DESCRIPTION:
!-----------------------------------------------------------------------
!  This routine receives message from cpl7 driver
!
!    The following fields are always received from the coupler:
! 
!    o  taux   -- zonal wind stress (taux)                 (W/m2   )
!    o  tauy   -- meridonal wind stress (tauy)             (W/m2   )
!    o  snow   -- water flux due to snow                   (kg/m2/s)
!    o  rain   -- water flux due to rain                   (kg/m2/s)
!    o  evap   -- evaporation flux                         (kg/m2/s)
!    o  meltw  -- snow melt flux                           (kg/m2/s)
!    o  salt   -- salt                                     (kg(salt)/m2/s)
!    o  swnet  -- net short-wave heat flux                 (W/m2   )
!    o  sen    -- sensible heat flux                       (W/m2   )
!    o  lwup   -- longwave radiation (up)                  (W/m2   )
!    o  lwdn   -- longwave radiation (down)                (W/m2   )
!    o  melth  -- heat flux from snow&ice melt             (W/m2   )
!    o  ifrac  -- ice fraction
!    o  rofl   -- river runoff flux                        (kg/m2/s)
!    o  rofi   -- ice runoff flux                          (kg/m2/s)
! 
!    The following fields are sometimes received from the coupler,
!      depending on model options:
! 
!    o  pslv   -- sea-level pressure                       (Pa)
!    o  duu10n -- 10m wind speed squared                   (m^2/s^2)
!    o  co2prog-- bottom atm level prognostic co2
!    o  co2diag-- bottom atm level diagnostic co2
! 
!-----------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

  !real(r8)           , intent(inout) :: x2o(:,:)
  !logical (log_kind) , intent(in)    :: ldiag_cpl
  !integer (POP_i4)   , intent(out)   :: errorCode  ! returned error code

  !TODO: update the types of following params
  double precision, intent(inout) :: x2o(:,:)
  logical, intent(in)             :: ldiag_cpl
  integer, intent(out)            :: errorCode  ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------










!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_import

!***********************************************************************
!BOP
! !IROUTINE: ocn_export_mct
! !INTERFACE:

 subroutine ocn_export(o2x, ldiag_cpl, errorCode)   

! !DESCRIPTION:
!  This routine calls the routines necessary to send MOM6 fields to
!  the CCSM cpl7 driver
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

  !real(r8)           , intent(inout) :: o2x(:,:)
  !logical (log_kind) , intent(in)    :: ldiag_cpl
  !integer (POP_i4)   , intent(out)   :: errorCode  ! returned error code

  !TODO: update the types of following params
  double precision, intent(inout) :: o2x(:,:)
  logical, intent(in)             :: ldiag_cpl
  integer, intent(out)            :: errorCode  ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------








!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_export

!***********************************************************************


end module ocn_import_export

