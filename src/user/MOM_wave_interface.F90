!> Initial conditions and forcing for the single column model (SCM) CVmix
!! test set.
module MOM_wave_interface

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, allocate_forcing_type
use MOM_grid, only : ocean_grid_type
use MOM_verticalgrid, only: verticalGrid_type
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/), get_time,&
                             time_type_to_real
use MOM_variables, only : thermo_var_ptrs, surface
use data_override_mod, only : data_override_init, data_override
implicit none ; private

#include <MOM_memory.h>

public MOM_wave_interface_init
public Import_Stokes_Drift


!> Container for wave related parameters
type, public:: wave_parameters_CS ;
private
  logical :: UseWaves    !< True to Compute Wave parameters
  integer :: WaveMethod  !< Options for various wave methods
  integer :: SpecMethod  !< Options for various wave spectra
  integer :: NumBands    !< Number of wavenumber bands to recieve
  real ALLOCABLE_, dimension(:) :: WaveNum_Cen !Wavenumber bands for read/coupled
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_), public :: &
       Us_x ! Stokes drift (zonal) 
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_), public :: &
       Us_y ! Stokes drift (meridional) 
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) ::&
       LangNum !Langmuir number
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,:), public :: &
       STKx0_Bands
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,:), public :: &
       STKy0_Bands  

  logical :: dataoverrideisinitialized
end type

! This include declares and sets the variable "version".
#include "version_variable.h"

character(len=40)  :: mod = "MOM_wave_interface" ! This module's name.

! Switches needed in import_stokes_drift
integer, parameter :: NO_SCHEME = 0, FROMFILE = 1, DATAOVERRIDE =2,&
                      PARAMETRIC = 3, TESTPROF = 99;
integer, parameter :: ELFOUHAILY = 1

! For Test Prof
Real :: TP_STKX0, TP_STKY0, TP_WVL

!-------------------------------------------------------------
CONTAINS
!
!> Initializes parameters related to MOM_wave_interface
subroutine MOM_wave_interface_init(G,GV,param_file, CS)
  type(ocean_grid_type),                  intent(in)  :: G !< Grid structure
  type(verticalGrid_type),                intent(in)  :: GV!< Vertical grid structure
  type(param_file_type),                  intent(in)  :: param_file !< Input parameter structure
  type(wave_parameters_CS),              pointer     :: CS
  ! Local variables
  integer :: is, ie, js, je, isd, ied, jsd, jed, isdB, iedB, jsdB, jedB, nz

  ! I/O
  character*(13) :: TMPSTRING1,TMPSTRING2
  character*(10), parameter :: NULL_STRING = "EMPTYEMPTY"
  character*(10), parameter :: PARAMETRIC_STRING = "PARAMETRIC" 
  character*(10), parameter :: FROMFILE_STRING = "FROM_FILE"
  character*(13), parameter :: DATAOVERRIDE_STRING = "DATA_OVERRIDE"
  character*(10), parameter :: TESTPROF_STRING = "TEST_PROF"
  character*(10), parameter :: ELF97_STRING = "ELF97"
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isdB = G%isdB ; iedB = G%iedB ; jsdB = G%jsdB ; jedB = G%jedB
  if (associated(CS)) then
     call MOM_error(WARNING, "wave_interface_init called with an associated"//&
                             "control structure.")
     return
  endif
  
  allocate(CS)

  call log_version(param_file, mod, version)
  call get_param(param_file,mod,"USE_WAVES",CS%UseWaves, &
                 'Main switch to use wave input', units='',default=.false.)

  if (CS%UseWaves) then 
     ! 1. Get Wave Method and write to integer WaveMethod
     call get_param(param_file,mod,"WAVE_METHOD",TMPSTRING1, &
                    'Choice of wave method, valid options...',units='',&
                    default=NULL_STRING)
     select case (TRIM(TMPSTRING1))
       case (NULL_STRING)! No Waves
          CS%WaveMethod = NO_SCHEME
          Print*,'You did not specify a wave method, so no waves are used.'
       case (FROMFILE_STRING)! From File
          CS%WaveMethod = FROMFILE
       case (DATAOVERRIDE_STRING)! From File
          CS%WaveMethod = DATAOVERRIDE
      case (PARAMETRIC_STRING)! Parameteric Spc
          CS%WaveMethod = PARAMETRIC
          call get_param(param_file,mod,"SPECTRUM_CHOICE",TMPSTRING2, &
                         'Choice of empirical wave spectrum, valid...',&
                         units='',default=NULL_STRING)
          select case(TRIM(TMPSTRING2))
            case (NULL_STRING)
               CS%WaveMethod = NO_SCHEME
               Print*,' '
               Print*,'You did not specify an empirical wave spectrum,'
               Print*,' so no waves are used.'
            case (ELF97_STRING)
               CS%SpecMethod = ELFOUHAILY
          endselect
        case (TESTPROF_STRING)
           CS%WaveMethod = TESTPROF
           call get_param(param_file,mod,"TP_STKX_SURF",TP_STKX0,&
                          'Surface Stokes (x) for test profile',&
                          units='m/s',default=0.1)
           call get_param(param_file,mod,"TP_STKY_SURF",TP_STKY0,&
                          'Surface Stokes (y) for test profile',&
                          units='m/s',default=0.0)
           call get_param(param_file,mod,"TP_WVL",TP_WVL,&
                          units='m',default=50.0)
     endselect
     ! 2. Allocate and initialize
     !    Stokes drift
     ALLOC_ (CS%Us_x(isdB:IedB,jsd:jed,nz)) ; CS%Us_x(:,:,:) = 0.0
     ALLOC_ (CS%Us_y(isd:Ied,jsdB:jedB,nz)) ; CS%Us_y(:,:,:) = 0.0
     !    Langmuir number
     ALLOC_ (CS%LangNum(isd:ied,jsd:jed)) ; CS%LangNum(:,:) = 1e10
  endif

  !/BGRTEMP{
  print*,' '
  print*,'-----------------------------------------------'
  print*,'You chose this wave method: ',CS%WaveMethod
  if(CS%WaveMethod==PARAMETRIC) then
     print*,'You chose this specrum: ',CS%SpecMethod
  endif
  if (CS%WaveMethod==TESTPROF) then
     print*,' '
     print*,'You chose the following for the test profile'
     print*,'--------------------------------------------'
     print*,'Surface Stk X [m/s]: ',TP_STKX0
     print*,'Surface Stk Y [m/s]: ',TP_STKY0
     print*,'Mean Wavelength [m]: ',TP_WVL
  endif
  print*,'-----------------------------------------------'
  print*,' '
  !\BGRTEMP}

  !----
  !/This will be moved, but for now call to construct Stokes drift profile here
  !----
  !Call Import_Stokes_Drift(G,GV,CS)

end subroutine MOM_wave_interface_init
!/
!/
!/
! Constructs the Stokes Drift profile on the model grid based on 
! desired coupling options
subroutine Import_Stokes_Drift(G,GV,Day,DT,CS,h)
  type(wave_parameters_CS),              pointer       :: CS
  type(ocean_grid_type),                  intent(in)   :: G !< Grid structure
  type(verticalGrid_type),                intent(in)   :: GV!< Vertical grid structure
  type(time_type), intent(in)                          :: Day
  type(time_type), intent(in)                          :: DT
 real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h
  ! local variables
  real    :: Top, MidPoint, Bottom
  real    :: DecayScale
  type(time_type) :: Day_Center
  integer :: ii, jj, kk
  integer :: is, ie, js, je, isd, ied, jsd, jed, isdB, iedB, jsdB, jedB, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isdB = G%isdB ; iedB = G%iedB ; jsdB = G%jsdB ; jedB = G%jedB

  Day_Center = Day + DT/2

  !print*,' '
  !print*,'Into Import_Stokes_Drift'

  if (CS%WaveMethod==TESTPROF) then
     DecayScale = 12.5663706/TP_WVL !4pi
     !print*,'Test Profile Construction'
     Bottom = 0.0
     MidPoint = 0.0
     !print*,'DecayScale:',TP_WVL,TP_STKX0,TP_STKY0
     do kk=1, nz
        Top = Bottom
        !****************************************************
        !NOTE THIS H WILL NOT BE CORRECT FOR NON-UNIFORM GRID
        MidPoint = Bottom - GV%H_to_m * h(1,1,kk)/2.
        Bottom = Bottom - GV%H_to_m * h(1,1,kk)
        do ii=isdB,iedB
           do jj=jsd,jed
              CS%Us_x(ii,jj,kk) = TP_STKX0 * EXP(MIDPOINT*DecayScale)
           enddo
        enddo
        do ii=isd,ied
           do jj=jsdB,jedB
              CS%Us_y(ii,jj,kk) = TP_STKY0 * EXP(MIDPOINT*DecayScale)
           enddo
        enddo
        !print*,MIDPOINT,CS%US_x(1,1,kk),CS%Us_y(1,1,kk)
     enddo     
  elseif (CS%WaveMethod==DATAOVERRIDE) then
     call Stokes_Drift_by_data_override(day_center,G,GV,CS,h)
  else!Keep this else, fallback to 0 Stokes drift
     do ii=isdB,iedB
           do jj=jsd,jed
              CS%Us_x(ii,jj,kk) = 0
           enddo
        enddo
        do ii=isd,ied
           do jj=jsdB,jedB
              CS%Us_y(ii,jj,kk) = 0
           enddo
        enddo
  endif

  !print*,'Out of Import_Stokes_Drift'
  !print*,' '
 
  !\BGRTEMP{
  !print*,' '
  !print*,'***********************************************'
  !print*,'** Made it to end of MOM_wave_interface_init **'
  !print*,'** Now stopping test...                      **'
  !print*,'***********************************************'
  !print*,' '
  !stop
  !\BGRTEMP}

end subroutine Import_Stokes_Drift
!
subroutine Stokes_Drift_by_data_override(day_center,G,GV,CS,h)
  use NETCDF
  type(time_type),             intent(in)  :: day_center
  type(wave_parameters_CS),    pointer     :: CS
  type(ocean_grid_type),       intent(in)  :: G !< Grid structure
  type(verticalGrid_type),     intent(in)  :: GV!< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h
  ! local variables
  real    :: Top, MidPoint, Bottom
  real    :: DecayScale
  integer :: i, j, is_in, ie_in, js_in, je_in

  integer, dimension(4) :: start, count, dims, dim_id 
  character(len=12)  :: dim_name(4)
  character(20) :: varname, filename, varread
  integer :: rcode, ncid, varid, id, ndims
  is_in = G%isc - G%isd + 1
  ie_in = G%iec - G%isd + 1
  js_in = G%jsc - G%jsd + 1
  je_in = G%jec - G%jsd + 1

  !/temp
  CS%dataOverrideIsInitialized = .false.

  if (.not.CS%dataOverrideIsInitialized) then

    call data_override_init(Ocean_domain_in=G%Domain%mpp_domain)
    CS%dataOverrideIsInitialized = .true.
    
    ! Read in number of wavenumber bands in file to set number to be read in
    filename = 'StkSpec1d.nc'
    varread = 'wavenumber'

    rcode = NF90_OPEN(trim(filename), NF90_NOWRITE, ncid)
    if (rcode .ne. 0) call MOM_error(FATAL,"error opening file "//trim(filename)//&
         " in MOM_wave_interface.")

    rcode = NF90_INQ_VARID(ncid, varread, varid)
    if (rcode .ne. 0) call MOM_error(FATAL,"error finding variable "//trim(varread)//&
         " in file "//trim(filename)//" in MOM_wave_interface.")

    rcode = NF90_INQUIRE_VARIABLE(ncid, varid, ndims=ndims, dimids=dims)
    if (rcode .ne. 0) call MOM_error(FATAL,'error inquiring dimensions MOM_wave_interface.')

    rcode = NF90_INQUIRE_DIMENSION(ncid, dims(1), dim_name(1), len=id)
    if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 1 data for "// &
         trim(varread)//" in file "// trim(filename)//" in MOM_wave_interface.")

    rcode = NF90_INQ_VARID(ncid, dim_name(1), dim_id(1))
    if (rcode .ne. 0) call MOM_error(FATAL,"error finding variable "//trim(dim_name(1))//&
         " in file "//trim(filename)//" in MOM_wave_interace.")

    ! Allocating size of wavenumber bins
    ALLOC_ ( CS%WaveNum_Cen(1:id) ) ; CS%WaveNum_Cen(:)=0.0


    ! Reading wavenumber bins
    start = 1; count = 1; count(1) = id
    rcode = NF90_GET_VAR(ncid, dim_id(1), CS%WaveNum_Cen, start, count)
  if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 1 values for var_name "// &
                trim(varread)//",dim_name "//trim(dim_name(1))//" in file "// trim(filename)//" in MOM_wave_interface")

    print*,CS%WaveNum_Cen

    print*,'******************'
    print*,'End of NetCDF Read'
 endif
  
  !do b=1,CS%NumBands
  !   write(varname,'(A15),(I0)') 'StokesX_Band',b
  !   call data_override('OCN',trim(varname), CS%STKx0_Bands(:,:,b), day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
  !   write(varname,'(A15),(I0)') 'StokesY_Band',b
  !   call data_override('OCN',trim(varname), CS%STKy0_Bands(:,:,b), day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
  !enddo

  print*,'End of Stokes Drift By Data Override.'
  stop

end subroutine Stokes_Drift_by_Data_Override
end module MOM_wave_interface
