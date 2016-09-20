!> This file contains the routines to read in and interpret wave data.
!> At present, the capabilities include reading in and returning Stokes drift.
module MOM_wave_interface

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : diag_ctrl
use MOM_domains,       only : pass_var, pass_vector, AGRID
use MOM_domains,       only : To_South, To_West, To_All
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, allocate_forcing_type
use MOM_grid, only : ocean_grid_type
use MOM_verticalgrid, only: verticalGrid_type
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/), get_time,&
                             time_type_to_real,real_to_time_type
use MOM_variables, only : thermo_var_ptrs, surface
use data_override_mod, only : data_override_init, data_override
implicit none ; private

#include <MOM_memory.h>

public MOM_wave_interface_init
public Import_Stokes_Drift
public StokesMixing
public CoriolisStokes
public Waves_end

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
  !^ Number of wavenumber/frequency partitions to receive
  !   - This needs to match the number of bands provided
  !     via the either coupling or file.
  !   - A check will be added to be sure it is correct 
  !     as this module nears completion.
  integer, public :: PartitionMode
  ! = 0 if partitions are described by wavenumber
  ! = 1 if partitions are described by frequency
  integer, public :: StkLevelMode=1
  ! = 0 if mid-point value of Stokes drift is used
  ! = 1 if average value of Stokes drift over level is used.
  real ALLOCABLE_, dimension(:), public :: &
   WaveNum_Cen,&   !< Wavenumber bands for read/coupled
   Freq_Cen        !< Frequency bands for read/coupled
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

  type(time_type), pointer, public :: Time ! A pointer to the ocean model's clock.
  type(diag_ctrl), pointer, public :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.

  integer, public :: id_StokesDrift_x, id_StokesDrift_y

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
subroutine MOM_wave_interface_init(time,G,GV,param_file, CS, diag )
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),                  intent(in)  :: G !< Grid structure
  type(verticalGrid_type),                intent(in)  :: GV!< Vertical grid structure
  type(param_file_type),                  intent(in)  :: param_file !< Input parameter structure
  type(wave_parameters_CS),              pointer     :: CS !< Wave parameter control structure
  type(diag_ctrl), target, intent(inout) :: diag
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

  CS%diag => diag
  CS%Time => Time

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

  return

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
  real    :: CMN_FAC, WN
  real, parameter :: PI=3.14159265359
  type(time_type) :: Day_Center
  integer :: ii, jj, kk, b, iim1, jjm1

  Day_Center = Day + DT/2

  !print*,' '
  !print*,'Into Import_Stokes_Drift'
  !print*,time_type_to_real(Day),time_type_to_real(Day_center)
  if (CS%WaveMethod==TESTPROF) then
     DecayScale = 4.*PI/TP_WVL !4pi
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
     do ii=G%isdB,G%iedB
        do jj=G%jsd,G%jed
           do b=1,CS%NumBands
              if (CS%PartitionMode==0) then
                 !In wavenumber we are averaging over (small) level
                 CMN_FAC = (1.0 - EXP(-0.01*2*CS%WaveNum_Cen(b))) / (0.01)/&
                      (2*CS%WaveNum_Cen(b))
              elseif (CS%PartitionMode==1) then
                 !In frequency we are not averaging over level and taking top
                 CMN_FAC = 1.0
              endif
              CS%US0_x(ii,jj)=CS%US0_x(ii,jj) + CS%STKx0(ii,jj,b) * CMN_FAC
           enddo

           bottom = 0.0
           do kk=1, G%ke
              Top = Bottom
              iim1 = max(ii-1,1)
              MidPoint = Bottom - GV%H_to_m * (h(ii,jj,kk)+h(iim1,jj,kk))/4.
              Bottom = Bottom - GV%H_to_m *  (h(ii,jj,kk)+h(iim1,jj,kk))/2.
              do b=1,CS%NumBands
                 if (CS%PartitionMode==0) then
                    !In wavenumber we are averaging over level
                    CMN_FAC =  (EXP(TOP*2.*CS%WaveNum_Cen(b))- &
                         EXP(BOTTOM*2.*CS%WaveNum_Cen(b))) / (Top-Bottom)/&
                         (2.*CS%WaveNum_Cen(b))
                 elseif (CS%PartitionMode==1) then
                    if (CS%StkLevelMode==0) then
                       !! Take the value at the midpoint
                       CMN_FAC = EXP(MidPoint*2.*(2.*PI*CS%Freq_Cen(b))**2/ &
                            GV%g_Earth)
                    elseif (CS%StkLevelMode==1) then
                       !! Use a numerical integration and then
                       !! divide by layer thickness
                       WN=(2.*PI*CS%Freq_Cen(b))**2
                       CMN_FAC = (exp(2.*WN*Top)-exp(2.*WN*Bottom)) &
                            /(2.*WN)/(Top-Bottom)
                    endif
                 endif
                 CS%US_x(ii,jj,kk)=CS%US_x(ii,jj,kk) + CS%STKx0(ii,jj,b) *&
                      CMN_FAC
              enddo
           enddo
        enddo
     enddo

     do ii=G%isd,G%ied
        do jj=G%jsdB,G%jedB
           do b=1,CS%NumBands
               if (CS%PartitionMode==0) then
                 !In wavenumber we are averaging over (small) level
                 CMN_FAC = (1.0 - EXP(-0.01*2*CS%WaveNum_Cen(b))) / (0.01)/&
                      (2*CS%WaveNum_Cen(b))
              elseif (CS%PartitionMode==1) then
                 !In frequency we are not averaging over level and taking top
                 CMN_FAC = 1.0
              endif
              CS%US0_y(ii,jj)=CS%US0_y(ii,jj) + CS%STKy0(ii,jj,b) * CMN_FAC
           enddo

           bottom = 0.0
           do kk=1, G%ke
              Top = Bottom
              jjm1 = max(jj-1,1)
              MidPoint = Bottom - GV%H_to_m * (h(ii,jj,kk)+h(ii,jjm1,kk))/4.
              Bottom = Bottom - GV%H_to_m *  (h(ii,jj,kk)+h(ii,jjm1,kk))/2.
              do b=1,CS%NumBands
                 if (CS%PartitionMode==0) then
                    !In wavenumber we are averaging over level
                    CMN_FAC =  (EXP(TOP*2.*CS%WaveNum_Cen(b))- &
                         EXP(BOTTOM*2.*CS%WaveNum_Cen(b))) / (Top-Bottom)/&
                         (2.*CS%WaveNum_Cen(b))
                 elseif (CS%PartitionMode==1) then
                    if (CS%StkLevelMode==0) then
                       !! Take the value at the midpoint
                       CMN_FAC = EXP(MidPoint*2.*(2.*PI*CS%Freq_Cen(b))**2/ &
                            GV%g_Earth)
                    elseif (CS%StkLevelMode==1) then
                       !! Use a numerical integration and then
                       !! divide by layer thickness
                       WN=(2.*PI*CS%Freq_Cen(b))**2
                       CMN_FAC = (exp(2.*WN*Top)-exp(2.*WN*Bottom)) &
                            /(2.*WN)/(Top-Bottom)
                    endif
                 endif
                 CS%US_y(ii,jj,kk)=CS%US_y(ii,jj,kk) + CS%STKy0(ii,jj,b) *&
                      CMN_FAC
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

!  print*,'USx'
!  print*,CS%US0_x(3,3),CS%Us_x(3,3,1)
!  print*,'USyx'
!  print*,CS%US0_x(3,3),CS%Us_x(3,3,1)


end subroutine Import_Stokes_Drift
!
subroutine Stokes_Drift_by_data_override(day_center,G,GV,CS)
  use NETCDF
  type(time_type),             intent(in)  :: day_center
  type(wave_parameters_CS),    pointer     :: CS
  type(ocean_grid_type),       intent(in)  :: G !< Grid structure
  type(verticalGrid_type),     intent(in)  :: GV!< Vertical grid structure
  ! local variables
  real    :: temp_x(SZI_(G),SZJ_(G)) ! Pseudo-zonal and psuedo-meridional
  real    :: temp_y(SZI_(G),SZJ_(G)) ! Stokex drift of band at h-points, in m/s
  real    :: Top, MidPoint, Bottom
  real    :: DecayScale
  integer :: b
  integer :: i, j

  integer, dimension(4) :: start, count, dims, dim_id 
  character(len=12)  :: dim_name(4)
  character(20) :: varname, filename, varread1, varread2
  integer :: rcode_fr, rcode_wn, ncid, varid_fr, varid_wn, id, ndims

  if (.not.CS%dataOverrideIsInitialized) then
    call data_override_init(Ocean_domain_in=G%Domain%mpp_domain)
    CS%dataOverrideIsInitialized = .true.

    ! Read in number of wavenumber bands in file to set number to be read in
    ! Hardcoded filename/variables
    filename = 'TMP.nc'
    varread1 = 'wavenumber' !Old method gives wavenumber
    varread2 = 'frequency'  !New method gives frequency
    rcode_wn = NF90_OPEN(trim(filename), NF90_NOWRITE, ncid)
    if (rcode_wn .ne. 0) then
       call MOM_error(FATAL,"error opening file "//trim(filename)//&
            " in MOM_wave_interface.")
    endif

    rcode_wn = NF90_INQ_VARID(ncid, varread1, varid_wn)
    rcode_fr = NF90_INQ_VARID(ncid, varread2, varid_fr)

    if (rcode_wn .ne. 0 .and. rcode_fr .ne. 0) then
       call MOM_error(FATAL,"error finding variable "//trim(varread1)//&
         " or "//trim(varread2)//" in file "//trim(filename)//" in MOM_wave_interface.")
    
    elseif (rcode_wn.eq.0) then
       ! wavenumbers found:
       CS%PartitionMode=0
       rcode_wn = NF90_INQUIRE_VARIABLE(ncid, varid_wn, ndims=ndims, &
            dimids=dims)
       if (rcode_wn .ne. 0) then
          call MOM_error(FATAL, &
               'error inquiring dimensions MOM_wave_interface.')
       endif
       rcode_wn = NF90_INQUIRE_DIMENSION(ncid, dims(1), dim_name(1), len=id)
       if (rcode_wn .ne. 0) then
          call MOM_error(FATAL,"error reading dimension 1 data for "// &
               trim(varread1)//" in file "// trim(filename)//          &
               " in MOM_wave_interface.")
       endif
       rcode_wn = NF90_INQ_VARID(ncid, dim_name(1), dim_id(1))
       if (rcode_wn .ne. 0) then
          call MOM_error(FATAL,"error finding variable "//trim(dim_name(1))//&
            " in file "//trim(filename)//" in MOM_wave_interace.")
       endif
       ! Allocating size of wavenumber bins
       ALLOC_ ( CS%WaveNum_Cen(1:id) ) ; CS%WaveNum_Cen(:)=0.0
    elseif (rcode_fr.eq.0) then
       ! frequencies found:
       CS%PartitionMode=1
       rcode_fr = NF90_INQUIRE_VARIABLE(ncid, varid_fr, ndims=ndims, &
            dimids=dims)
       if (rcode_fr .ne. 0) then
          call MOM_error(FATAL,&
               'error inquiring dimensions MOM_wave_interface.')
       endif
       rcode_fr = NF90_INQUIRE_DIMENSION(ncid, dims(1), dim_name(1), len=id)
       if (rcode_fr .ne. 0) then
          call MOM_error(FATAL,"error reading dimension 1 data for "// &
               trim(varread2)//" in file "// trim(filename)// &
               " in MOM_wave_interface.")
       endif
       rcode_fr = NF90_INQ_VARID(ncid, dim_name(1), dim_id(1))
       if (rcode_fr .ne. 0) then
          call MOM_error(FATAL,"error finding variable "//trim(dim_name(1))//&
               " in file "//trim(filename)//" in MOM_wave_interace.")
       endif
       ! Allocating size of frequency bins
       ALLOC_ ( CS%Freq_Cen(1:id) ) ; CS%Freq_Cen(:)=0.0
    endif


    ! Allocating size of wavenumber bins
    
    ALLOC_ ( CS%STKx0(G%isdB:G%iedB,G%jsd:G%jed,1:id)) ; CS%STKx0(:,:,:) = 0.0
    ALLOC_ ( CS%STKy0(G%isd:G%ied,G%jsdB:G%jedB,1:id)) ; CS%STKy0(:,:,:) = 0.0
    

    ! Reading wavenumber bins/Frequencies
    start = 1; count = 1; count(1) = id
    if (CS%PartitionMode==0) then
       rcode_wn = NF90_GET_VAR(ncid, dim_id(1), CS%WaveNum_Cen, start, count)
       if (rcode_wn .ne. 0) then
          call MOM_error(FATAL,&
               "error reading dimension 1 values for var_name "// &
               trim(varread1)//",dim_name "//trim(dim_name(1))//  &
               " in file "// trim(filename)//" in MOM_wave_interface")
       endif
       CS%NUMBANDS = ID
    elseif (CS%PartitionMode==1) then
       rcode_fr = NF90_GET_VAR(ncid, dim_id(1), CS%Freq_Cen, start, count)
       if (rcode_fr .ne. 0) then
          call MOM_error(FATAL,&
               "error reading dimension 1 values for var_name "// &
               trim(varread2)//",dim_name "//trim(dim_name(1))//  &
               " in file "// trim(filename)//" in MOM_wave_interface")
       endif
       CS%NUMBANDS = ID
    endif

  endif

  do b=1,CS%NumBands
    temp_x(:,:)=0.0;temp_y(:,:)=0.0;
    varname = '                    '
    write(varname,"(A3,I0)")'Usx',b
    call data_override('OCN',trim(varname), temp_x, day_center)
    varname = '                    '
    write(varname,'(A3,I0)')'Usy',b
    call data_override('OCN',trim(varname), temp_y, day_center)
    ! Disperse into halo on h-grid
    call pass_vector(temp_x, temp_y, G%Domain, To_All, AGRID)
    !Filter land values
    do j = G%jsd,G%jed ; do I = G%Isd,G%Ied
       if (abs(temp_x(i,j)).gt.10. .or. abs(temp_y(i,j)).gt.10. ) then
          ! Assume land-mask and zero out
          temp_x(i,j)=0.0
          temp_y(i,j)=0.0
       endif
    enddo; enddo
    ! Interpolate to u/v grids
    do j = G%jsc,G%jec ; do I = G%IscB,G%IecB
       CS%STKx0(I,j,b) = 0.5 * (temp_x(i,j) + temp_x(i+1,j))
    enddo; enddo
    do j = G%JscB,G%JecB ; do i = G%isc,G%iec
       CS%STKy0(i,J,b) = 0.5 * (temp_y(i,j) + temp_y(i,j+1))
    enddo; enddo
    ! Disperse into halo on u/v grids
    call pass_vector(CS%STKx0(:,:,b),CS%STKy0(:,:,b), G%Domain, To_ALL)
    !print*,'maxXandY: ',maxval(abs(temp_x)),maxval(abs(temp_y))
    !print*,'maxXandY: ',maxval(abs(CS%STKx0(:,:,b))),maxval(abs(CS%STKy0(:,:,b)))
  enddo

return

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

!> Clear pointers, deallocate memory
subroutine Waves_end(CS)
!/
  type(wave_parameters_CS), pointer :: CS !< Control structure
!/
  if (allocated(CS%WaveNum_Cen)) then; DEALLOC_( CS%WaveNum_Cen ); endif
  if (allocated(CS%Freq_Cen))    DEALLOC_( CS%Freq_Cen )
  if (allocated(CS%Us_x))        DEALLOC_( CS%Us_x )
  if (allocated(CS%Us_y))        DEALLOC_( CS%Us_y )
  if (allocated(CS%LangNum))     DEALLOC_( CS%LangNum )
  if (allocated(CS%STKx0))       DEALLOC_( CS%STKx0 )
  if (allocated(CS%STKy0))       DEALLOC_( CS%STKy0 )
  if (allocated(CS%KvS))         DEALLOC_( CS%KvS )
  if (allocated(CS%Us0_y))       DEALLOC_( CS%Us0_y )
  if (allocated(CS%Us0_x))       DEALLOC_( CS%Us0_x )
!/
  DEALLOC_( CS )
!/
  return
end subroutine Waves_end
end module MOM_wave_interface
