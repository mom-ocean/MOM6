!> This file contains the routines to read in and interpret wave data.
!> At present, the capabilities include reading in and returning Stokes drift.
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
public StokesMixing
public CoriolisStokes

!> Container for wave related parameters
type, public:: wave_parameters_CS ;
private
  logical :: UseWaves
  !^ True to enable this module
  logical, public :: LagrangianMixing
  !^ True if Stokes drift is present and mixing
  !  should be applied to Lagrangian current
  !  (mean current + Stokes drift) 
  logical, public :: StokesMixing
  !^ True if Stokes drift is present and mixing
  !  should be applied directly to Stokes current
  !  (with separate mixing parameter for Eulerian
  !  mixing contribution)
  logical, public :: CoriolisStokes
  !^ True if Stokes drift is present and Coriolis
  !  Stokes acceleration of mean current 
  !  should be applied.
  logical, public :: StokesInKappaShear=.false.!This doesn't work well.
  !^ True if Stokes drift is present and 
  !  should be added into Kappa Shear mixing.
  logical, public :: LangmuirEnhanceW   
  !^ True if turbulent velocity scales should be
  !  enhanced due to Langmuir mixing.
  logical, public :: LangmuirEnhanceVt2 
  !^ True if unresolved turbulent velocity scale
  ! should be enhanced due to Langmuir mixing.
  logical, public :: LangmuirEnhanceK   
  !^ True if diffusivity and viscosity should be
  !  enhanced due to presence of Langmuir mixing.
  logical, public :: StokesShearInRIb   
  !^ True if Stokes drift profile should be included
  !  in current shear calculation for bulk Richardson number.
  logical, public ::SurfaceStokesInRIb  
  !^ True if surface Stokes drift should be used
  !  to enhance the denominator of bulk Richardson number.
  integer :: WaveMethod                 
  !^ Options for various wave methods
  !    Valid (tested) choices are:
  !     - TEST_PROF
  !     - DATA_OVERRIDE
  !    In prep include:
  !     - Parametric (See SpecMethod)
  integer :: SpecMethod  
  !^ Options for various Parametric wave spectra
  !     IN PREP
  integer, public :: NumBands    
  !^ Number of wavenumber bands to receive
  !   - This needs to match the number of bands provided
  !     via the either coupling or file.
  !   - A check will be added to be sure it is correct 
  !     as this module nears completion.
  real ALLOCABLE_, dimension(:), public :: &
   WaveNum_Cen   !< Wavenumber bands for read/coupled
  real ALLOCABLE_, dimension( NIMEMB_, NJMEM_,NKMEM_), public :: &
       Us_x !< Stokes drift profile (zonal) 
  real ALLOCABLE_, dimension( NIMEM_, NJMEMB_,NKMEM_), public :: &
       Us_y !< Stokes drift profile (meridional) 
  real ALLOCABLE_, dimension( NIMEM_, NJMEM_), public ::         &
       LangNum, & !< Langmuir number (directionality factored later)
       US0_x,   & !< Surface Stokes Drift in x
       US0_y      !< Surface Stokes Drift in y
  real ALLOCABLE_, dimension( NIMEMB_, NJMEM_,NKMEM_), public :: &
       STKx0 !< Stokes Drift spectrum in x
  real ALLOCABLE_, dimension( NIMEM_, NJMEMB_,NKMEM_), public :: &
       STKy0 !< Stokes Drift spectrum in y 
  real ALLOCABLE_, dimension( NIMEM_, NJMEM_,NKMEM_), public :: &
       KvS !< Viscosity for Stokes Drift shear 
  logical :: dataoverrideisinitialized
end type wave_parameters_CS

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
  type(wave_parameters_CS),              pointer     :: CS !< Wave parameter control structure
  ! Local variables

  ! I/O
  character*(13) :: TMPSTRING1,TMPSTRING2
  character*(10), parameter :: NULL_STRING = "EMPTYEMPTY"
  character*(10), parameter :: PARAMETRIC_STRING = "PARAMETRIC" 
  character*(10), parameter :: FROMFILE_STRING = "FROM_FILE"
  character*(13), parameter :: DATAOVERRIDE_STRING = "DATA_OVERRIDE"
  character*(10), parameter :: TESTPROF_STRING = "TEST_PROF"
  character*(10), parameter :: ELF97_STRING = "ELF97"

  if (associated(CS)) then
     call MOM_error(WARNING, "wave_interface_init called with an associated"//&
                             "control structure.")
     return
  endif
  
  allocate(CS)

  ! Add any initializations needed here
  CS%dataOverrideIsInitialized = .false.
  ! Here we assume the only way to get here is with UseWaves enabled.
  ! Therefore, we do not read from input file again.
  CS%UseWaves=.true.

  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "LAGRANGIAN_MIXING", CS%LagrangianMixing, &
       "Flag to use Lagrangian Mixing of momentum", units="", &
       Default=.false.)
  call get_param(param_file, mod, "STOKES_MIXING", CS%StokesMixing, &
       "Flag to use Stokes Mixing of momentum", units="", &
       Default=.false.)  
  call get_param(param_file, mod, "CORIOLIS_STOKES", CS%CoriolisStokes, &
       "Flag to use Coriolis Stokes acceleration", units="", &
       Default=.false.)  
  call get_param(param_file, mod, "LANGMUIR_ENHANCE_W", CS%LangmuirEnhanceW, &
       'Flag for Langmuir turbulence enhancement of turbulent'//&
       'velocity scale.', units="", Default=.false.) 
  call get_param(param_file, mod, "LANGMUIR_ENHANCE_VT2", CS%LangmuirEnhanceVt2, &
       'Flag for Langmuir turbulence enhancement of Vt2 in KPP'//&
       'bulk Richardson Number.', units="", Default=.false.) 
  call get_param(param_file, mod, "LANGMUIR_ENHANCE_K", CS%LangmuirEnhanceK, &
       'Flag for Langmuir turbulence enhancement of turbulent'//&
       'mixing coefficients.', units="", Default=.false.) 
 call get_param(param_file, mod, "STOKES_IN_RIB", CS%StokesShearInRIb, &
       'Flag for using Stokes drift profile in RIb calculation.'&
       , units="", Default=.false.) 
 call get_param(param_file, mod, "SURFACE_STOKES_IN_RIB", CS%SurfaceStokesInRIb, &
       'Flag for using surface Stokes drift in RIb calculation.'&
       , units="", Default=.false.) 
  if ( (CS%LagrangianMixing.or.CS%LangmuirEnhanceW) .and. (.not.CS%UseWaves)) then
     call MOM_error(FATAL,"MOM_vert_friction(visc): "// &
          "LagrangianMixing and WaveEnhancedDiff cannot"//&
          "be called without USE_WAVES = .true.")
  endif
  ! Return if not using waves.  Should not have gotten here to begin with.
  if (.not.CS%UseWaves) return

  ! 1. Get Wave Method and write to integer WaveMethod
  call get_param(param_file,mod,"WAVE_METHOD",TMPSTRING1,   &
       "Choice of wave method, valid options include: \n"// &
       "  DATA_OVERRIDE - Get wave info from NetCDF \n"//   &
       "  TEST_PROF - Monochromatic spectrum\n"//           &
       "     Define: TP_STKX_SURF,TP_STKy_SURF, & TP_WVL.", &
       units='', default=NULL_STRING)
  select case (TRIM(TMPSTRING1))
  case (NULL_STRING)! No Waves
     CS%WaveMethod = NO_SCHEME
     Print*,'You did not specify a wave method, so no waves are used.'
  case (FROMFILE_STRING)! From File
     CS%WaveMethod = FROMFILE
     PRINT*,'Read from file not Coded, try DATA_OVERRIDE'
  case (DATAOVERRIDE_STRING)! From DataOverride NetCDF
     CS%WaveMethod = DATAOVERRIDE
  case (PARAMETRIC_STRING)! Use Parameteric Spectrum
     CS%WaveMethod = PARAMETRIC
     PRINT*,'Parametric spectrum not coded'
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
  ALLOC_ (CS%Us_x(G%isdB:G%IedB,G%jsd:G%jed,G%ke)) ; CS%Us_x(:,:,:) = 0.0
  ALLOC_ (CS%Us_y(G%isd:G%Ied,G%jsdB:G%jedB,G%ke)) ; CS%Us_y(:,:,:) = 0.0
  !    Langmuir number
  ALLOC_ (CS%LangNum(G%isc:G%iec,G%jsc:G%jec)) ; CS%LangNum(:,:) = 1e10
  ALLOC_ (CS%US0_x(G%isdB:G%iedB,G%jsd:G%jed)) ; CS%US0_x(:,:) = 0.
  ALLOC_ (CS%US0_y(G%isd:G%ied,G%jsdB:G%jedB)) ; CS%US0_y(:,:) = 0.  
  ! Viscosity for Stokes drift
  ALLOC_ (CS%KvS(G%isd:G%Ied,G%jsd:G%jed,G%ke)) ; CS%KvS(:,:,:) = 1.e-6

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

end subroutine MOM_wave_interface_init
!/
!/
!/
! Constructs the Stokes Drift profile on the model grid based on 
! desired coupling options
subroutine Import_Stokes_Drift(G,GV,Day,DT,CS,h,FLUXES)
  type(wave_parameters_CS),              pointer        :: CS
  type(ocean_grid_type),                  intent(in)    :: G !< Grid structure
  type(verticalGrid_type),                intent(in)    :: GV!< Vertical grid structure
  type(time_type), intent(in)                           :: Day
  type(time_type), intent(in)                           :: DT
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h
  type(forcing), intent(in)                             :: FLUXES
  ! local variables
  real    :: USy20pct, USx20pct, H20pct
  real    :: Top, MidPoint, Bottom
  real    :: DecayScale
  type(time_type) :: Day_Center
  integer :: ii, jj, kk, b, iim1, jjm1

  Day_Center = Day + DT/2

  !print*,' '
  !print*,'Into Import_Stokes_Drift'

  if (CS%WaveMethod==TESTPROF) then
     DecayScale = 12.5663706/TP_WVL !4pi
     do ii=G%isdB,G%iedB
        do jj=G%jsd,G%jed
           iim1=max(1,ii-1)
           Bottom = 0.0
           MidPoint = 0.0
           do kk=1, G%ke
              Top = Bottom
              MidPoint = Bottom - GV%H_to_m * (h(ii,jj,kk)+h(iim1,jj,kk))/4.
              Bottom = Bottom - GV%H_to_m * (h(ii,jj,kk)+h(iim1,jj,kk))/2.
              CS%Us_x(ii,jj,kk) = TP_STKX0 * EXP(MIDPOINT*DecayScale)
           enddo
        enddo
     enddo
     do ii=G%isd,G%ied
        do jj=G%jsdB,G%jedB
           jjm1=max(1,jj-1)
           Bottom = 0.0
           MidPoint = 0.0
           do kk=1, G%ke
              Top = Bottom
              MidPoint = Bottom - GV%H_to_m * (h(ii,jj,kk)+h(ii,jjm1,kk))/4.
              Bottom = Bottom - GV%H_to_m * (h(ii,jj,kk)+h(ii,jjm1,kk))/2.
              CS%Us_y(ii,jj,kk) = TP_STKY0 * EXP(MIDPOINT*DecayScale)
           enddo
        enddo
     enddo
  elseif (CS%WaveMethod==DATAOVERRIDE) then
       call Stokes_Drift_by_data_override(day_center,G,GV,CS)
       CS%Us_x(:,:,:) = 0.0
       CS%Us_y(:,:,:) = 0.0
       CS%Us0_x(:,:) = 0.0
       CS%Us0_y(:,:) = 0.0
     ! ---------------------------------------------------------|
     ! This computes the average Stokes drift based on the      |
     !  analytical integral over the layer divided by the layer |
     !  thickness.                                              |
     ! ---------------------------------------------------------|
     do ii=G%isdB,G%iedB
        do jj=G%jsd,G%jed
           do b=1,CS%NumBands
              CS%US0_x(ii,jj)=CS%US0_x(ii,jj) + CS%STKx0(ii,jj,b) *&
                   (1.0 - EXP(-0.01*2*CS%WaveNum_Cen(b))) / (0.01)/&
                   (2*CS%WaveNum_Cen(b))
           enddo
           bottom = 0.0
           do kk=1, G%ke
              Top = Bottom
              iim1 = max(ii-1,1)
              MidPoint = Bottom - GV%H_to_m * (h(ii,jj,kk)+h(iim1,jj,kk))/4.
              Bottom = Bottom - GV%H_to_m *  (h(ii,jj,kk)+h(iim1,jj,kk))/2.
              do b=1,CS%NumBands
                 CS%US_x(ii,jj,kk)=CS%US_x(ii,jj,kk) + CS%STKx0(ii,jj,b) *&
                      (EXP(TOP*2*CS%WaveNum_Cen(b))- &
                      EXP(BOTTOM*2*CS%WaveNum_Cen(b))) / (Top-Bottom)/&
                       (2*CS%WaveNum_Cen(b))
              enddo
           enddo
        enddo
     enddo

     do ii=G%isd,G%ied
        do jj=G%jsdB,G%jedB
           do b=1,CS%NumBands
              CS%US0_y(ii,jj)=CS%US0_y(ii,jj) + CS%STKy0(ii,jj,b) *&
                   (1.0 - EXP(-0.01*2*CS%WaveNum_Cen(b))) / (0.01)/&
                   (2*CS%WaveNum_Cen(b))
           enddo
           bottom = 0.0
           do kk=1, G%ke
              Top = Bottom
              jjm1 = max(jj-1,1)
              MidPoint = Bottom - GV%H_to_m * (h(ii,jj,kk)+h(ii,jjm1,kk))/4.
              Bottom = Bottom - GV%H_to_m *  (h(ii,jj,kk)+h(ii,jjm1,kk))/2.
              do b=1,CS%NumBands
                 CS%US_y(ii,jj,kk)=CS%US_y(ii,jj,kk) + CS%STKy0(ii,jj,b) *&
                      (EXP(TOP*2*CS%WaveNum_Cen(b))- &
                      EXP(BOTTOM*2*CS%WaveNum_Cen(b))) / (Top-Bottom)/&
                       (2*CS%WaveNum_Cen(b))
              enddo
           enddo
        enddo
     enddo

  else!Keep this else, fallback to 0 Stokes drift
     do ii=G%isdB,G%iedB
           do jj=G%jsd,G%jed
              CS%Us_x(ii,jj,kk) = 0.
           enddo
        enddo
        do ii=G%isd,G%ied
           do jj=G%jsdB,G%jedB
              CS%Us_y(ii,jj,kk) = 0.
           enddo
        enddo
  endif

end subroutine Import_Stokes_Drift
!
subroutine Stokes_Drift_by_data_override(day_center,G,GV,CS)
  use NETCDF
  type(time_type),             intent(in)  :: day_center
  type(wave_parameters_CS),    pointer     :: CS
  type(ocean_grid_type),       intent(in)  :: G !< Grid structure
  type(verticalGrid_type),     intent(in)  :: GV!< Vertical grid structure
  ! local variables
  real    :: Top, MidPoint, Bottom
  real    :: DecayScale
  integer :: b
  integer :: i, j

  integer, dimension(4) :: start, count, dims, dim_id 
  character(len=12)  :: dim_name(4)
  character(20) :: varname, filename, varread
  integer :: rcode, ncid, varid, id, ndims

  if (.not.CS%dataOverrideIsInitialized) then
    print*,'into init'
    call data_override_init(Ocean_domain_in=G%Domain%mpp_domain)
    print*,'out of init'
    CS%dataOverrideIsInitialized = .true.
    
    ! Read in number of wavenumber bands in file to set number to be read in
    ! Hardcoded filename/variables
    filename = 'StkSpec.nc'
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
    ALLOC_ ( CS%STKx0(G%isdB:G%iedB,G%jsd:G%jed,1:id)) ; CS%STKx0(:,:,:) = 0.0
    ALLOC_ ( CS%STKy0(G%isd:G%ied,G%jsdB:G%jedB,1:id)) ; CS%STKy0(:,:,:) = 0.0
    

    ! Reading wavenumber bins
    start = 1; count = 1; count(1) = id
    rcode = NF90_GET_VAR(ncid, dim_id(1), CS%WaveNum_Cen, start, count)
    if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 1 values for var_name "// &
         trim(varread)//",dim_name "//trim(dim_name(1))//" in file "// trim(filename)//" in MOM_wave_interface")

    CS%NUMBANDS = ID
  endif
  
  do b=1,CS%NumBands
     varname = '                    '
     write(varname,"(A3,I0)")'Usx',b
     call data_override('OCN',trim(varname), CS%STKx0(:,:,b), day_center)
     varname = '                    '
     write(varname,'(A3,I0)')'Usy',b
     call data_override('OCN',trim(varname), CS%STKy0(:,:,b), day_center)     
  enddo
  
  !Brandon: Hacking to update HALO until properly resolved
  do i=G%isdB,G%iedB
     do j=G%isd,G%ied
        if (i.lt.G%iscB) then
           CS%STKX0(i,j,:)=CS%STKX0(G%iscB,j,:)
        elseif (i.gt.G%iecB) then
           CS%STKX0(i,j,:)=CS%STKX0(G%iecB,j,:)
        endif
     enddo
  enddo
  do i=G%isd,G%ied
     do j=G%isdB,G%iedB
        if (j.lt.G%jscB) then
           CS%STKY0(i,j,:)=CS%STKY0(i,G%iscB,:)
        elseif (j.gt.G%jecB) then
           CS%STKY0(i,j,:)=CS%STKY0(i,G%jecB,:)
        endif
     enddo
  enddo
  
end subroutine Stokes_Drift_by_Data_Override
!/
!/
!/
subroutine StokesMixing(G, GV, DT, h, u, v, WAVES, FLUXES)
  ! Arguments
  type(ocean_grid_type),                  intent(in)    :: G              !< Ocean grid
  type(verticalGrid_type),                intent(in)    :: GV             !< Ocean vertical grid
  real, intent(in)                                         :: Dt             !< Time step of MOM6 [s] for GOTM turbulence solver
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h              !< Layer/level thicknesses (units of H)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: u              !< Velocity i-component (m/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)    :: v              !< Velocity j-component (m/s)
  type(Wave_parameters_CS), pointer                         :: Waves  !< Surface wave related control structure.
  type(forcing), intent(in)                             :: FLUXES
  ! Local variables
  REAL :: dTauUp, dTauDn, DVel
  INTEGER :: i,j,k

! This is a very poor way to do Stokes mixing.
!  Cannot separate Stokes/Eulerian mixing due to boundary condition.
!  This is really just a temporary attempt...

  do k = 1, G%ke
     do j = G%jscB, G%jecB
        do i = G%iscB, G%iecB
           if (k.eq.1) then
              dTauUp = 0.!FLUXES%taux(i,j)/1000.!Convert????
              dTauDn =  0.5*(WAVES%Kvs(i,j,k+1)+WAVES%Kvs(i+1,j,k+1))*&
                   (waves%us_x(i,j,k)-waves%us_x(i,j,k+1))&
                   /(GV%H_to_m *0.5*(h(i,j,k)+h(i,j,k+1)) )
           elseif (k.lt.G%ke-1) then
              dTauUp =   0.5*(waves%Kvs(i,j,k)+waves%Kvs(i+1,j,k))*&
                   (waves%us_x(i,j,k-1)-waves%us_x(i,j,k))&
                   /(GV%H_to_m *0.5*(h(i,j,k-1)+h(i,j,k)) )
              dTauDn =  0.5*(waves%Kvs(i,j,k+1)+waves%Kvs(i+1,j,k+1))*&
                   (waves%us_x(i,j,k)-waves%us_x(i,j,k+1))&
                   /(GV%H_to_m *0.5*(h(i,j,k)+h(i,j,k+1)) )
           elseif (k.eq.G%ke) then
              dTauUp =   0.5*(waves%Kvs(i,j,k)+waves%Kvs(i+1,j,k))*&
                   (waves%us_x(i,j,k-1)-waves%us_x(i,j,k))&
                   /(GV%H_to_m *0.5*(h(i,j,k-1)+h(i,j,k)) )
              dTauDn = 0.0!FLUXES%taux
           endif
           DVel = (dTauUp-dTauDn) / (GV%H_to_m *h(i,j,k)) * DT
           !if (i.eq.3 .and. j.eq.3 .and. k.eq.1) then
           !   print*,u(i,j,k),(DVel),dtaudn
           !endif
           u(i,j,k) = u(i,j,k)+DVel
        enddo
     enddo
  enddo
     

  do k = 1, G%ke
     do j = G%jscB, G%jecB
        do i = G%iscB, G%iecB
           if (k.eq.1) then
              dTauUp = 0.!FLUXES%tauy(i,j)/1000.!Convert????
              dTauDn = 0.5*(waves%Kvs(i,j,k+1)+waves%Kvs(i,j+1,k+1))&
                   *(waves%us_y(i,j,k)-waves%us_y(i,j,k+1))&
                   /(GV%H_to_m *0.5*(h(i,j,k)+h(i,j,k+1)) )
           elseif (k.lt.G%ke-1) then
              dTauUp =   0.5*(waves%Kvs(i,j,k)+waves%Kvs(i,j+1,k))*&
                   (waves%us_y(i,j,k-1)-waves%us_y(i,j,k))&
                   /(GV%H_to_m *0.5*(h(i,j,k-1)+h(i,j,k)) )
              dTauDn =  0.5*(waves%Kvs(i,j,k+1)+waves%Kvs(i,j+1,k+1))*&
                   (waves%us_y(i,j,k)-waves%us_y(i,j,k+1))&
                   /(GV%H_to_m *0.5*(h(i,j,k)+h(i,j,k+1)) )
           elseif (k.eq.G%ke) then
              dTauUp =   0.5*(waves%Kvs(i,j,k)+waves%Kvs(i,j+1,k))*&
                   (waves%us_y(i,j,k-1)-waves%us_y(i,j,k))&
                   /(GV%H_to_m *0.5*(h(i,j,k-1)+h(i,j,k)) )
              dTauDn = 0.0!FLUXES%tauy
           endif
           DVel = (dTauUp-dTauDn) / (GV%H_to_m *h(i,j,k)) * DT
           !if (i.eq.3 .and. j.eq.3 .and. k.eq.1) then
           !   print*,v(i,j,k),(DVel),dtaudn
           !endif
           v(i,j,k) = v(i,j,k)+DVel
        enddo
     enddo
  enddo

end subroutine StokesMixing

subroutine CoriolisStokes(G, GV, DT, h, u, v, WAVES)
  ! Arguments
  type(ocean_grid_type),                  intent(in)    :: G              !< Ocean grid
  type(verticalGrid_type),                intent(in)    :: GV             !< Ocean vertical grid
  real, intent(in)                                         :: Dt             !< Time step of MOM6 [s] for GOTM turbulence solver
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h              !< Layer/level thicknesses (units of H)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: u              !< Velocity i-component (m/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)    :: v              !< Velocity j-component (m/s)
  type(Wave_parameters_CS), pointer                         :: Waves  !< Surface wave related control structure.

  ! Local variables
  REAL :: DVel
  INTEGER :: i,j,k

  do k = 1, G%ke
     do j = G%jscB, G%jecB
        do i = G%iscB, G%iecB
           DVel = 0.25*(WAVES%us_y(i,j+1,k)+WAVES%us_y(i-1,j+1,k))*G%CoriolisBu(i,j+1) +  0.25*(WAVES%us_y(i,j,k)+WAVES%us_y(i-1,j,k))*G%CoriolisBu(i,j)
           u(i,j,k) = u(i,j,k)+DVEL
        enddo
     enddo
  enddo
     
  do k = 1, G%ke
     do j = G%jscB, G%jecB
        do i = G%iscB, G%iecB
           DVel = 0.25*(WAVES%us_x(i+1,j,k)+WAVES%us_x(i+1,j-1,k))*G%CoriolisBu(i+1,j) +  0.25*(WAVES%us_x(i,j,k)+WAVES%us_x(i,j-1,k))*G%CoriolisBu(i,j)
           v(i,j,k) = v(i,j,k)-DVEL
        enddo
     enddo
  enddo
end subroutine CoriolisStokes
end module MOM_wave_interface
