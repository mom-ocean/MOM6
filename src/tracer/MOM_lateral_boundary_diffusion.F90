!> Calculates and applies diffusive fluxes as a parameterization of lateral mixing (non-neutral) by
!! mesoscale eddies near the top and bottom (to be implemented) boundary layers of the ocean.

module MOM_lateral_boundary_diffusion

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,             only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_domains,               only : pass_var
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_diag_mediator,         only : post_data, register_diag_field
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type, log_param
use MOM_file_parser,           only : openParameterBlock, closeParameterBlock
use MOM_io,                    only : file_exists, field_size, MOM_read_data, slasher, field_exists
use MOM_grid,                  only : ocean_grid_type
use MOM_remapping,             only : remapping_CS, initialize_remapping
use MOM_remapping,             only : extract_member_remapping_CS, build_reconstructions_1d
use MOM_remapping,             only : average_value_ppoly, remappingSchemesDoc, remappingDefaultScheme
use MOM_remapping,             only : remapping_core_h
use MOM_regridding,            only : check_grid_def
use MOM_tracer_registry,       only : tracer_registry_type, tracer_type
use MOM_unit_scaling,          only : unit_scale_type
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_CVMix_KPP,             only : KPP_get_BLD, KPP_CS
use MOM_energetic_PBL,         only : energetic_PBL_get_MLD, energetic_PBL_CS
use MOM_diabatic_driver,       only : diabatic_CS, extract_diabatic_member
use MOM_string_functions,      only : extract_integer, extract_real, extractWord
use iso_fortran_env,           only : stdout=>output_unit, stderr=>error_unit

implicit none ; private

public near_boundary_unit_tests, lateral_boundary_diffusion, lateral_boundary_diffusion_init
public boundary_k_range

! Private parameters to avoid doing string comparisons for bottom or top boundary layer
integer, public, parameter :: SURFACE = -1 !< Set a value that corresponds to the surface bopundary
integer, public, parameter :: BOTTOM  = 1  !< Set a value that corresponds to the bottom boundary
#include <MOM_memory.h>

!> Sets parameters for lateral boundary mixing module.
type, public :: lbd_CS ; private
  integer :: deg                             !< Degree of polynomial reconstruction
  integer :: nk                              !< Number of layers in dz_top
  integer :: surface_boundary_scheme         !< Which boundary layer scheme to use
                                             !! 1. ePBL; 2. KPP
  logical :: linear                          !< If True, apply a linear transition at the base/top of the boundary.
                                             !! The flux will be fully applied at k=k_min and zero at k=k_max.
  real, dimension(:), allocatable  :: dz_top !< top vertical grid to remap the state before applying lateral diffusion
  real, dimension(:), allocatable  :: dz_bot !< bot vertical grid to remap the state before applying lateral diffusion
  type(remapping_CS)              :: remap_CS                     !< Control structure to hold remapping configuration
  type(KPP_CS),           pointer :: KPP_CSp => NULL()            !< KPP control structure needed to get BLD
  type(energetic_PBL_CS), pointer :: energetic_PBL_CSp => NULL()  !< ePBL control structure needed to get BLD
  type(diag_ctrl), pointer :: diag => NULL()                      !< A structure that is used to
                                                                  !! regulate the timing of diagnostic output.
end type lbd_CS

! This include declares and sets the variable "version".
#include "version_variable.h"
character(len=40) :: mdl = "MOM_lateral_boundary_diffusion"       !< Name of this module

contains

!> Initialization routine that reads runtime parameters and sets up pointers to other control structures that might be
!! needed for lateral boundary diffusion.
logical function lateral_boundary_diffusion_init(Time, G, param_file, diag, diabatic_CSp, CS)
  type(time_type), target,          intent(in)    :: Time          !< Time structure
  type(ocean_grid_type),            intent(in)    :: G             !< Grid structure
  type(param_file_type),            intent(in)    :: param_file    !< Parameter file structure
  type(diag_ctrl), target,          intent(inout) :: diag          !< Diagnostics control structure
  type(diabatic_CS),                pointer       :: diabatic_CSp  !< KPP control structure needed to get BLD
  type(lbd_CS),                     pointer       :: CS            !< Lateral boundary mixing control structure

  ! local variables
  character(len=80)  :: string, varName    ! Temporary strings
  character(len=200) :: inputdir, fileName ! Temporary strings
  character(len=320) :: message            ! Temporary strings
  character(len=12)  :: expected_units     ! Temporary strings
  integer :: ke, nk ! Number of levels in the LBD and native grids, respectively
  logical :: boundary_extrap ! controls if boundary extrapolation is used in the LBD code
  logical :: ierr
  real    :: tmpReal
  integer :: nzf(4)
  real, dimension(:), allocatable :: z_max  ! Maximum interface depths [H ~> m or kg m-2] or other
                                            ! units depending on the coordinate
  if (ASSOCIATED(CS)) then
    call MOM_error(FATAL, "lateral_boundary_diffusion_init called with associated control structure.")
    return
  endif

  ! Log this module and master switch for turning it on/off
  call log_version(param_file, mdl, version, &
       "This module implements lateral diffusion of tracers near boundaries")
  call get_param(param_file, mdl, "USE_LATERAL_BOUNDARY_DIFFUSION", lateral_boundary_diffusion_init, &
                 "If true, enables the lateral boundary tracer's diffusion module.", &
                 default=.false.)

  if (.not. lateral_boundary_diffusion_init) then
    return
  endif

  allocate(CS)
  CS%diag => diag
  call extract_diabatic_member(diabatic_CSp, KPP_CSp=CS%KPP_CSp)
  call extract_diabatic_member(diabatic_CSp, energetic_PBL_CSp=CS%energetic_PBL_CSp)

  CS%surface_boundary_scheme = -1
  !GMM, uncomment below
!  if ( .not. ASSOCIATED(CS%energetic_PBL_CSp) .and. .not. ASSOCIATED(CS%KPP_CSp) ) then
!    call MOM_error(FATAL,"Lateral boundary diffusion is true, but no valid boundary layer scheme was found")
!  endif

  ! Read all relevant parameters and write them to the model log.
  call get_param(param_file, mdl, "LBD_LINEAR_TRANSITION", CS%linear, &
                 "If True, apply a linear transition at the base/top of the boundary. \n"//&
                 "The flux will be fully applied at k=k_min and zero at k=k_max.", default=.false.)
  call get_param(param_file, mdl, "LBD_BOUNDARY_EXTRAP", boundary_extrap, &
                 "Use boundary extrapolation in LBD code", &
                 default=.false.)
  call get_param(param_file, mdl, "LBD_REMAPPING_SCHEME", string, &
                 "This sets the reconstruction scheme used "//&
                 "for vertical remapping for all variables. "//&
                 "It can be one of the following schemes: "//&
                 trim(remappingSchemesDoc), default=remappingDefaultScheme)
  call initialize_remapping( CS%remap_CS, string, boundary_extrapolation = boundary_extrap )
  call extract_member_remapping_CS(CS%remap_CS, degree=CS%deg)
  ! set dz_top
  call get_param(param_file, mdl, "LBD_DIAG_COORD_TOP", string, &
                 "Determines how to specify the vertical resolution "//&
                 "to apply lateral diffusion near the surface. Valid options are:\n"//&
                 " PARAM       - use the vector-parameter LBD_DZ_TOP \n"//&
                 " UNIFORM[:N] - uniformly distributed\n"//&
                 " FILE:string - read from a file. The string specifies\n"//&
                 "               the filename and variable name, separated\n"//&
                 "               by a comma or space, e.g. FILE:lev.nc,dz\n"//&
                 "               or FILE:lev.nc,interfaces=zw\n",&
                 default="UNIFORM:500,500")
  message = "The distribution of vertical resolution used to \n"//&
            "apply lateral diffusion near boundaries."
  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".", do_not_log=.true.)
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "NK", nk, &
                 "The number of model layers.", units="nondim", fail_if_missing=.true., &
                 do_not_log=.true.)
  if (index(trim(string),'UNIFORM')==1) then
    call get_param(param_file, "MOM", "MAXIMUM_DEPTH", tmpReal, &
                 "The maximum depth of the ocean.", units="m", default=4000.0, do_not_log=.true.)
    if (len_trim(string)==7) then
      ke = nk ! Use model nk by default
    elseif (index(trim(string),'UNIFORM:')==1 .and. len_trim(string)>8) then
      ! Format is "UNIFORM:N" or "UNIFORM:N,MAX_DEPTH"
      ke = extract_integer(string(9:len_trim(string)),'',1)
      tmpReal = extract_real(string(9:len_trim(string)),',',2,missing_value=tmpReal)
    else
      call MOM_error(FATAL,trim(mdl)//', lateral_boundary_diffusion_init: '// &
          'Unable to interpret "'//trim(string)//'".')
    endif
    allocate(CS%dz_top(ke))
    CS%dz_top(:) = tmpReal / real(ke)
    call log_param(param_file, mdl, "!LBD_DZ_TOP", CS%dz_top, &
                   trim(message), units='m')
  elseif (trim(string)=='PARAM') then
    ke = nk ! Use model nk by default
    allocate(CS%dz_top(ke))
    call get_param(param_file, mdl, 'LBD_DZ_TOP', CS%dz_top, &
                   trim(message), units='m', fail_if_missing=.true.)
  elseif (index(trim(string),'FILE:')==1) then
    ! FILE:filename,var_name is assumed to be reading level thickness variables
    ! FILE:filename,interfaces=var_name reads positions
    if (string(6:6)=='.' .or. string(6:6)=='/') then
      ! If we specified "FILE:./xyz" or "FILE:/xyz" then we have a relative or absolute path
      fileName = trim( extractWord(trim(string(6:80)), 1) )
    else
      ! Otherwise assume we should look for the file in INPUTDIR
      fileName = trim(inputdir) // trim( extractWord(trim(string(6:80)), 1) )
    endif
    if (.not. file_exists(fileName)) call MOM_error(FATAL,trim(mdl)//", lateral_boundary_diffusion_init: "// &
            "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")

    varName = trim( extractWord(trim(string(6:)), 2) )
    if (len_trim(varName)==0) then
      if (field_exists(fileName,'dz')) then; varName = 'dz'
      else ;  call MOM_error(FATAL,trim(mdl)//", lateral_boundary_diffusion_init: "// &
                    "Coordinate variable (dz) not specified and none could be guessed.")
      endif
    endif
    expected_units = 'meters'
    if (index(trim(varName),'interfaces=')==1) then
      varName=trim(varName(12:))
      call check_grid_def(filename, varName, expected_units, message, ierr)
      if (ierr) call MOM_error(FATAL,trim(mdl)//", lateral_boundary_diffusion_init: "//&
                  "Unsupported format in grid definition '"//trim(filename)//"'. Error message "//trim(message))
      call field_size(trim(fileName), trim(varName), nzf)
      ke = nzf(1)-1
      allocate(CS%dz_top(ke))
      allocate(z_max(ke+1))
      call MOM_read_data(trim(fileName), trim(varName), z_max)
      CS%dz_top(:) = abs(z_max(1:ke) - z_max(2:ke+1))
      deallocate(z_max)
    else
      ! Assume reading resolution
      call field_size(trim(fileName), trim(varName), nzf)
      ke = nzf(1)
      allocate(CS%dz_top(ke))
      call MOM_read_data(trim(fileName), trim(varName), CS%dz_top)
    endif
    call log_param(param_file, mdl, "!LBD_DZ_TOP", CS%dz_top, &
                   trim(message), units='m')
  else
    call MOM_error(FATAL,trim(mdl)//", lateral_boundary_diffusion_init: "// &
      "Unrecognized coordinate configuration"//trim(string))
  endif
  CS%nk = ke
  ! TODO: set dz_bot
  CS%dz_bot(:) = 1.0
end function lateral_boundary_diffusion_init

!> Driver routine for calculating lateral diffusive fluxes near the top and bottom boundaries.
!! Two different methods are available:
!! Method 1: more straight forward, diffusion is applied layer by layer using only information
!! from neighboring cells.
!! Method 2: lower order representation, calculate fluxes from bulk layer integrated quantities.
subroutine lateral_boundary_diffusion(G, GV, US, h, Coef_x, Coef_y, dt, Reg, CS)
  type(ocean_grid_type),                intent(inout) :: G   !< Grid type
  type(verticalGrid_type),              intent(in)    :: GV  !< ocean vertical grid structure
  type(unit_scale_type),                intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                                        intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)),    intent(in)    :: Coef_x !< dt * Kh * dy / dx at u-points [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G)),    intent(in)    :: Coef_y !< dt * Kh * dx / dy at v-points [L2 ~> m2]
  real,                                 intent(in)    :: dt     !< Tracer time step * I_numitts
                                                                !! (I_numitts in tracer_hordiff)
  type(tracer_registry_type),           pointer       :: Reg    !< Tracer registry
  type(lbd_CS),                         pointer       :: CS     !< Control structure for this module

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: hbl                           !< bnd. layer depth [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G),CS%deg+1) :: ppoly0_coefs !< Coefficients of polynomial
  real, dimension(SZI_(G),SZJ_(G),SZK_(G),2)        :: ppoly0_E     !< Edge values from reconstructions
  real, dimension(SZK_(G),CS%deg+1)                 :: ppoly_S      !< Slopes from reconstruction (placeholder)
  !real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: uFlx       !< Zonal flux of tracer [conc m^3]
  real, dimension(SZIB_(G),SZJ_(G),CS%nk)   :: uFlx        !< Zonal flux of tracer in z-space [conc m^3]
  real, dimension(SZIB_(G),SZJ_(G))         :: uFLx_bulk   !< Total calculated bulk-layer u-flux for the tracer
  !real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: vFlx       !< Meridional flux of tracer [conc m^3]
  real, dimension(SZI_(G),SZJB_(G),CS%nk)   :: vFlx        !< Meridional flux of tracer in z-space [conc m^3]
  real, dimension(SZI_(G),SZJB_(G))         :: vFlx_bulk   !< Total calculated bulk-layer v-flux for the tracer
  real, dimension(SZIB_(G),SZJ_(G))         :: uwork_2d    !< Layer summed u-flux transport
  real, dimension(SZI_(G),SZJB_(G))         :: vwork_2d    !< Layer summed v-flux transport
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: tendency    !< tendency array for diag in the zgrid
  real, dimension(SZI_(G),SZJ_(G),CS%nk)    :: tracer_z    !< Tracer in the zgrid
  real, dimension(SZI_(G),SZJ_(G))          :: tendency_2d !< depth integrated content tendency for diagn
  type(tracer_type), pointer                :: tracer => NULL() !< Pointer to the current tracer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: tracer_old  !< local copy of the initial tracer concentration,
                                                           !! only used to compute tendencies.
  real, dimension(SZI_(G),SZJ_(G),CS%nk)    :: tracer_z_old!< Copy of the initial tracer concentration in z-space

  integer :: remap_method !< Reconstruction method
  integer :: i,j,k,m      !< indices to loop over
  real    :: Idt          !< inverse of the time step [s-1]

  Idt = 1./dt
  hbl(:,:) = 100.
  hbl(4:6,:) = 100.
  !if (ASSOCIATED(CS%KPP_CSp)) call KPP_get_BLD(CS%KPP_CSp, hbl, G)
  !if (ASSOCIATED(CS%energetic_PBL_CSp)) call energetic_PBL_get_MLD(CS%energetic_PBL_CSp, hbl, G, US)

  call pass_var(hbl,G%Domain)
  do m = 1,Reg%ntr
    tracer => Reg%tr(m)
    tracer_z(:,:,:) = 0.0
    tracer_z_old(:,:,:) = 0.0
    ! for diagnostics
    if (tracer%id_lbdxy_conc > 0 .or. tracer%id_lbdxy_cont > 0 .or. tracer%id_lbdxy_cont_2d > 0) then
      tendency(:,:,:) = 0.0
      tracer_old(:,:,:) = 0.0
      ! copy initial tracer state so that the tendency can be computed
      tracer_old(:,:,:) = tracer%t(:,:,:)
    endif

    ! remap tracer to zgrid
    do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
      call remapping_core_h(CS%remap_cs, G%ke, h(i,j,:), tracer%t(i,j,:), CS%nk, CS%dz_top(:), tracer_z(i,j,:))
      !call build_reconstructions_1d( CS%remap_CS, G%ke, h(i,j,:), tracer%t(i,j,:), ppoly0_coefs(i,j,:,:), &
      !                               ppoly0_E(i,j,:,:), ppoly_S, remap_method, GV%H_subroundoff, GV%H_subroundoff)
    enddo ; enddo

    ! Diffusive fluxes in the i- and j-direction
    uFlx(:,:,:) = 0. ! z-space
    vFlx(:,:,:) = 0. ! z-space
    uFlx_bulk(:,:) = 0.
    vFlx_bulk(:,:) = 0.

    ! LBD layer by layer
    do j=G%jsc,G%jec
      do i=G%isc-1,G%iec
        if (G%mask2dCu(I,j)>0.) then
          !call fluxes_layer_method(SURFACE, GV%ke, CS%nk, CS%deg, h(I,j,:), h(I+1,j,:), hbl(I,j), hbl(I+1,j),  &
          !  G%areaT(I,j), G%areaT(I+1,j), tracer%t(I,j,:), tracer%t(I+1,j,:), ppoly0_coefs(I,j,:,:),    &
          !  ppoly0_coefs(I+1,j,:,:), ppoly0_E(I,j,:,:), ppoly0_E(I+1,j,:,:), remap_method, Coef_x(I,j), &
          !  uFlx(I,j,:), CS)
          call fluxes_layer_method1(SURFACE, CS%nk, hbl(I,j), hbl(I+1,j),  &
            G%areaT(I,j), G%areaT(I+1,j), tracer_z(I,j,:), tracer_z(I+1,j,:), &
            remap_method, Coef_x(I,j), uFlx(I,j,:), CS)
        endif
      enddo
    enddo
    do J=G%jsc-1,G%jec
      do i=G%isc,G%iec
        if (G%mask2dCv(i,J)>0.) then
          !call fluxes_layer_method(SURFACE, GV%ke, CS%nk, CS%deg, h(i,J,:), h(i,J+1,:), hbl(i,J), hbl(i,J+1),  &
          !  G%areaT(i,J), G%areaT(i,J+1), tracer%t(i,J,:), tracer%t(i,J+1,:), ppoly0_coefs(i,J,:,:),    &
          !  ppoly0_coefs(i,J+1,:,:), ppoly0_E(i,J,:,:), ppoly0_E(i,J+1,:,:), remap_method, Coef_y(i,J), &
          !  vFlx(i,J,:), CS)
          call fluxes_layer_method1(SURFACE, CS%nk, hbl(i,J), hbl(i,J+1),  &
            G%areaT(i,J), G%areaT(i,J+1), tracer_z(i,J,:), tracer_z(i,J+1,:), &
            remap_method, Coef_y(i,J), vFlx(i,J,:), CS)
        endif
      enddo
    enddo

    ! Update the tracer fluxes
    do k=1,CS%nk ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (G%mask2dT(i,j)>0.) then
        !tracer%t(i,j,k) = tracer%t(i,j,k) + (( (uFlx(I-1,j,k)-uFlx(I,j,k)) ) + ( (vFlx(i,J-1,k)-vFlx(i,J,k) ) ))* &
        !                  (G%IareaT(i,j)/( h(i,j,k) + GV%H_subroundoff))
        tracer_z(i,j,k) = tracer_z(i,j,k) + (( (uFlx(I-1,j,k)-uFlx(I,j,k)) ) + ( (vFlx(i,J-1,k)-vFlx(i,J,k) ) ))* &
                          (G%IareaT(i,j)/( CS%dz_top(k) + GV%H_subroundoff))
        ! difference between before/after diffusion in the zgrid
        tendency_z(i,j,k) = tracer_z(i,j,k) - tracer_z_old(i,j,k)
      endif
    enddo ; enddo ; enddo

    ! remap tracer "change" back to native grid
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      tracer_1d(:) = 0.0
      call remapping_core_h(CS%remap_cs, CS%nk, CS%dz_top, tendency_z(i,j,:), G%ke, h(i,j,:), tracer_1d(:))
      tracer%t(i,j,:) = tracer%t(i,j,:) + tracer_1d(:)
    enddo ; enddo

    if (tracer%id_lbdxy_conc > 0  .or. tracer%id_lbdxy_cont > 0 .or. tracer%id_lbdxy_cont_2d > 0 ) then
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
        tendency(i,j,k) = (tracer%t(i,j,k)-tracer_old(i,j,k)) * Idt
      enddo ; enddo ; enddo
    endif

    ! Post the tracer diagnostics
    !if (tracer%id_lbd_dfx>0)      call post_data(tracer%id_lbd_dfx, uFlx*Idt, CS%diag)
    !if (tracer%id_lbd_dfy>0)      call post_data(tracer%id_lbd_dfy, vFlx*Idt, CS%diag)
    if (tracer%id_lbd_dfx_2d>0) then
      uwork_2d(:,:) = 0.
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
        uwork_2d(I,j) = uwork_2d(I,j) + (uFlx(I,j,k) * Idt)
      enddo ; enddo ; enddo
      call post_data(tracer%id_lbd_dfx_2d, uwork_2d, CS%diag)
    endif

    if (tracer%id_lbd_dfy_2d>0) then
      vwork_2d(:,:) = 0.
      do k=1,GV%ke ; do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
        vwork_2d(i,J) = vwork_2d(i,J) + (vFlx(i,J,k) * Idt)
      enddo ; enddo ; enddo
      call post_data(tracer%id_lbd_dfy_2d, vwork_2d, CS%diag)
    endif

    ! post tendency of tracer content
    if (tracer%id_lbdxy_cont > 0) then
      call post_data(tracer%id_lbdxy_cont, tendency, CS%diag)
    endif

    ! post depth summed tendency for tracer content
    if (tracer%id_lbdxy_cont_2d > 0) then
      tendency_2d(:,:) = 0.
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        do k=1,GV%ke
          tendency_2d(i,j) = tendency_2d(i,j) + tendency(i,j,k)
        enddo
      enddo ; enddo
      call post_data(tracer%id_lbdxy_cont_2d, tendency_2d, CS%diag)
    endif

    ! post tendency of tracer concentration; this step must be
    ! done after posting tracer content tendency, since we alter
    ! the tendency array and its units.
    if (tracer%id_lbdxy_conc > 0) then
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
        tendency(i,j,k) =  tendency(i,j,k) / ( h(i,j,k) + GV%H_subroundoff )
      enddo ; enddo ; enddo
      call post_data(tracer%id_lbdxy_conc, tendency, CS%diag)
    endif

  enddo

end subroutine lateral_boundary_diffusion

!> Calculate the harmonic mean of two quantities
!! See \ref section_harmonic_mean.
real function harmonic_mean(h1,h2)
  real :: h1 !< Scalar quantity
  real :: h2 !< Scalar quantity
  if (h1 + h2 == 0.) then
    harmonic_mean = 0.
  else
    harmonic_mean = 2.*(h1*h2)/(h1+h2)
  endif
end function harmonic_mean

!> Find the k-index range corresponding to the layers that are within the boundary-layer region
subroutine boundary_k_range(boundary, nk, h, hbl, k_top, zeta_top, k_bot, zeta_bot)
  integer,             intent(in   ) :: boundary !< SURFACE or BOTTOM                       [nondim]
  integer,             intent(in   ) :: nk       !< Number of layers                        [nondim]
  real, dimension(nk), intent(in   ) :: h        !< Layer thicknesses of the column         [H ~> m or kg m-2]
  real,                intent(in   ) :: hbl      !< Thickness of the boundary layer         [H ~> m or kg m-2]
                                                 !! If surface, with respect to zbl_ref = 0.
                                                 !! If bottom, with respect to zbl_ref = SUM(h)
  integer,             intent(  out) :: k_top    !< Index of the first layer within the boundary
  real,                intent(  out) :: zeta_top !< Distance from the top of a layer to the intersection of the
                                                 !! top extent of the boundary layer (0 at top, 1 at bottom)  [nondim]
  integer,             intent(  out) :: k_bot    !< Index of the last layer within the boundary
  real,                intent(  out) :: zeta_bot !< Distance of the lower layer to the boundary layer depth
                                                 !! (0 at top, 1 at bottom)  [nondim]
  ! Local variables
  real :: htot ! Summed thickness [H ~> m or kg m-2]
  integer :: k
  ! Surface boundary layer
  if ( boundary == SURFACE ) then
    k_top = 1
    zeta_top = 0.
    htot = 0.
    k_bot = 1
    zeta_bot = 0.
    if (hbl == 0.) return
    if (hbl >= SUM(h(:))) then
      k_bot = nk
      zeta_bot = 1.
      return
    endif
    do k=1,nk
      htot = htot + h(k)
      if ( htot >= hbl) then
        k_bot = k
        zeta_bot = 1 - (htot - hbl)/h(k)
        return
      endif
    enddo
  ! Bottom boundary layer
  elseif ( boundary == BOTTOM ) then
    k_top = nk
    zeta_top = 1.
    k_bot = nk
    zeta_bot = 0.
    htot = 0.
    if (hbl == 0.) return
    if (hbl >= SUM(h(:))) then
      k_top = 1
      zeta_top = 1.
      return
    endif
    do k=nk,1,-1
      htot = htot + h(k)
      if (htot >= hbl) then
        k_top = k
        zeta_top = 1 - (htot - hbl)/h(k)
        return
      endif
    enddo
  else
    call MOM_error(FATAL,"Houston, we've had a problem in boundary_k_range")
  endif

end subroutine boundary_k_range

!> Calculate the lateral boundary diffusive fluxes using the layer by layer method.
!! See \ref section_method1
subroutine fluxes_layer_method1(boundary, nk, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, &
                              method, khtr_u, F_layer, CS)

  integer,                      intent(in   )    :: boundary !< Which boundary layer SURFACE or BOTTOM  [nondim]
  integer,                      intent(in   )    :: nk       !< Number of layers in the local z-grid    [nondim]
  real,                         intent(in   )    :: hbl_L    !< Thickness of the boundary boundary
                                                                       !! layer (left)              [H ~> m or kg m-2]
  real,                         intent(in   )    :: hbl_R    !< Thickness of the boundary boundary
                                                             !! layer (right)                       [H ~> m or kg m-2]
  real,                         intent(in   )    :: area_L   !< Area of the horizontal grid (left)  [L2 ~> m2]
  real,                         intent(in   )    :: area_R   !< Area of the horizontal grid (right) [L2 ~> m2]
  real, dimension(nk),        intent(in   )      :: phi_L    !< Tracer values (left)                [conc]
  real, dimension(nk),        intent(in   )      :: phi_R    !< Tracer values (right)               [conc]
  integer,                   intent(in   )       :: method   !< Method of polynomial integration    [nondim]
  real,                      intent(in   )       :: khtr_u   !< Horizontal diffusivities times delta t
                                                             !! at a velocity point [L2 ~> m2]
  real, dimension(nk),     intent(  out)         :: F_layer  !< Layerwise diffusive flux at U- or V-point in the local
                                                             !! z-grid [H L2 conc ~> m3 conc]
  type(lbd_CS),              pointer             :: CS       !< Lateral diffusion control structure
                                                             !! the boundary layer
  ! Local variables
  real                :: khtr_avg             !< Thickness-weighted diffusivity at the u-point     [m^2 s^-1]
                                              !! This is just to remind developers that khtr_avg should be
                                              !! computed once khtr is 3D.
  real                :: heff                 !< Harmonic mean of layer thicknesses           [H ~> m or kg m-2]
  real                :: inv_heff             !< Inverse of the harmonic mean of layer thicknesses
                                              !!  [H-1 ~> m-1 or m2 kg-1]
  real                :: phi_L_avg, phi_R_avg !< Bulk, thickness-weighted tracer averages (left and right column)
                                              !!                                                   [conc m^-3 ]
  real    :: htot                      !< Total column thickness [H ~> m or kg m-2]
  !real    :: heff_tot                  !< Total effective column thickness in the transition layer [m]
  integer :: k, k_bot_min, k_top_max   !< k-indices, min and max for bottom and top, respectively
  integer :: k_bot_max, k_top_min      !< k-indices, max and min for bottom and top, respectively
  integer :: k_bot_diff, k_top_diff    !< different between left and right k-indices for bottom and top, respectively
  integer :: k_top_L, k_bot_L          !< k-indices left native grid
  integer :: k_top_R, k_bot_R          !< k-indices right native grid
  real    :: zeta_top_L, zeta_top_R    !< distance from the top of a layer to the boundary
                                       !! layer depth in the native grid                  [nondim]
  real    :: zeta_bot_L, zeta_bot_R    !< distance from the bottom of a layer to the boundary
                                       !!layer depth in the native grid                   [nondim]
  real    :: h_work_L, h_work_R  !< dummy variables
  real    :: hbl_min             !< minimum BLD (left and right)                          [m]
  real    :: wgt                 !< weight to be used in the linear transition to the interior [nondim]
  real    :: a                   !< coefficient to be used in the linear transition to the interior [nondim]

  F_layer(:) = 0.0
  if (hbl_L == 0. .or. hbl_R == 0.) then
    return
  endif

  ! Calculate vertical indices containing the boundary layer in dz_top
  call boundary_k_range(boundary, nk, CS%dz_top, hbl_L, k_top_L, zeta_top_L, k_bot_L, zeta_bot_L)
  call boundary_k_range(boundary, nk, CS%dz_top, hbl_R, k_top_R, zeta_top_R, k_bot_R, zeta_bot_R)

  if (boundary == SURFACE) then
    k_bot_min = MIN(k_bot_L, k_bot_R)
    k_bot_max = MAX(k_bot_L, k_bot_R)
    k_bot_diff = (k_bot_max - k_bot_min)

    ! make sure left and right k indices span same range
    if (k_bot_min .ne. k_bot_L) then
      k_bot_L = k_bot_min
      zeta_bot_L = 1.0
    endif
    if (k_bot_min .ne. k_bot_R) then
      k_bot_R= k_bot_min
      zeta_bot_R = 1.0
    endif

    h_work_L = (CS%dz_top(k_bot_L) * zeta_bot_L)
    h_work_R = (CS%dz_top(k_bot_R) * zeta_bot_R)

    ! GMM, the following needs to be modified. We need to calculate ppoly0_E_L and ppoly0_coefs_L here...
    !phi_L_avg = average_value_ppoly( nk, phi_L_local, ppoly0_E_L, ppoly0_coefs_L, method, k_bot_L, 0., zeta_bot_L)
    !phi_R_avg = average_value_ppoly( nk, phi_R_local, ppoly0_E_R, ppoly0_coefs_R, method, k_bot_R, 0., zeta_bot_R)
    !heff = harmonic_mean(h_work_L, h_work_R)

    ! tracer flux where the minimum BLD intersets layer
    ! GMM, khtr_avg should be computed once khtr is 3D
    if ((CS%linear) .and. (k_bot_diff .gt. 1)) then
      ! apply linear decay at the base of hbl
      do k = k_bot_min-1,1,-1
        !heff = harmonic_mean(h_L(k), h_R(k))
        F_layer(k) = -(CS%dz_top(k) * khtr_u) * (phi_R(k) - phi_L(k))
      enddo
      htot = 0.0
      do k = k_bot_min+1,k_bot_max, 1
        htot = htot + CS%dz_top(k)
      enddo

      a = -1.0/htot
      htot = 0.0
      do k = k_bot_min,k_bot_max, 1
        !heff = harmonic_mean(h_L(k), h_R(k))
        wgt = (a*(htot + (CS%dz_top(k) * 0.5))) + 1.0
        F_layer(k) = -(CS%dz_top(k) * khtr_u) * (phi_R(k) - phi_L(k)) * wgt
        htot = htot + CS%dz_top(k)
      enddo
    else
      !!F_layer(k_bot_min) = -(heff * khtr_u) * (phi_R_avg - phi_L_avg)
      do k = k_bot_min-1,1,-1
        !heff = harmonic_mean(h_L(k), h_R(k))
        F_layer(k) = -(CS%dz_top(k) * khtr_u) * (phi_R(k) - phi_L(k))
      enddo
    endif
  endif

!  if (boundary == BOTTOM) then
!    ! TODO: GMM add option to apply linear decay
!    k_top_max = MAX(k_top_L, k_top_R)
!    ! make sure left and right k indices span same range
!    if (k_top_max .ne. k_top_L) then
!      k_top_L = k_top_max
!      zeta_top_L = 1.0
!    endif
!    if (k_top_max .ne. k_top_R) then
!      k_top_R= k_top_max
!      zeta_top_R = 1.0
!    endif
!
!    h_work_L = (CS%dz_bot(k_top_L) * zeta_top_L)
!    h_work_R = (CS%dz_bot(k_top_R) * zeta_top_R)
!
!    phi_L_avg = average_value_ppoly( nk, phi_L, ppoly0_E_L, ppoly0_coefs_L, method, k_top_L, 1.0-zeta_top_L, 1.0)
!    phi_R_avg = average_value_ppoly( nk, phi_R, ppoly0_E_R, ppoly0_coefs_R, method, k_top_R, 1.0-zeta_top_R, 1.0)
!    heff = harmonic_mean(h_work_L, h_work_R)
!
!    ! tracer flux where the minimum BLD intersets layer
!    F_layer(k_top_max) = (-heff * khtr_u) * (phi_R_avg - phi_L_avg)
!
!    do k = k_top_max+1,nk
!      heff = harmonic_mean(h_L(k), h_R(k))
!      F_layer(k) = -(heff * khtr_u) * (phi_R(k) - phi_L(k))
!    enddo
!  endif

end subroutine fluxes_layer_method1

!> Calculate the lateral boundary diffusive fluxes using the layer by layer method.
!! See \ref section_method1
subroutine fluxes_layer_method(boundary, nk, nk_z, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, &
                              ppoly0_coefs_L, ppoly0_coefs_R, ppoly0_E_L, ppoly0_E_R, method, khtr_u, &
                              F_layer, CS)

  integer,                      intent(in   )    :: boundary !< Which boundary layer SURFACE or BOTTOM  [nondim]
  integer,                      intent(in   )    :: nk       !< Number of layers in the native grid     [nondim]
  integer,                      intent(in   )    :: nk_z     !< Number of layers in the local z-grid    [nondim]
  integer,                      intent(in   )    :: deg      !< order of the polynomial reconstruction  [nondim]
  real, dimension(nk),          intent(in   )    :: h_L      !< Layer thickness (left)              [H ~> m or kg m-2]
  real, dimension(nk),          intent(in   )    :: h_R      !< Layer thickness (right)             [H ~> m or kg m-2]
  real,                         intent(in   )    :: hbl_L    !< Thickness of the boundary boundary
                                                                       !! layer (left)              [H ~> m or kg m-2]
  real,                         intent(in   )    :: hbl_R    !< Thickness of the boundary boundary
                                                             !! layer (right)                       [H ~> m or kg m-2]
  real,                         intent(in   )    :: area_L   !< Area of the horizontal grid (left)  [L2 ~> m2]
  real,                         intent(in   )    :: area_R   !< Area of the horizontal grid (right) [L2 ~> m2]
  real, dimension(nk),          intent(in   )    :: phi_L    !< Tracer values (left)                [conc]
  real, dimension(nk),          intent(in   )    :: phi_R    !< Tracer values (right)               [conc]
  real, dimension(nk,deg+1), intent(in   )    :: ppoly0_coefs_L !< Tracer reconstruction (left)  [conc]
  real, dimension(nk,deg+1), intent(in   )    :: ppoly0_coefs_R !< Tracer reconstruction (right) [conc]
  real, dimension(nk,2),     intent(in   )       :: ppoly0_E_L !< Polynomial edge values (left)     [nondim]
  real, dimension(nk,2),     intent(in   )       :: ppoly0_E_R !< Polynomial edge values (right)    [nondim]
  integer,                   intent(in   )       :: method   !< Method of polynomial integration    [nondim]
  real,                      intent(in   )       :: khtr_u   !< Horizontal diffusivities times delta t
                                                             !! at a velocity point [L2 ~> m2]
  real, dimension(nk_z),     intent(  out)       :: F_layer  !< Layerwise diffusive flux at U- or V-point in the local
                                                             !! z-grid [H L2 conc ~> m3 conc]
  type(lbd_CS),              pointer             :: CS       !< Lateral diffusion control structure
                                                             !! the boundary layer
  ! Local variables
  real, dimension(nk_z) :: phi_L_local        !< Tracer values (left) in the zgrid                 [conc]
  real, dimension(nk_z) :: phi_R_local        !< Tracer values (right) in the zgrid                [conc]
  real, dimension(nk) :: h_means              !< Calculate the layer-wise harmonic means           [H ~> m or kg m-2]
  real                :: khtr_avg             !< Thickness-weighted diffusivity at the u-point     [m^2 s^-1]
                                              !! This is just to remind developers that khtr_avg should be
                                              !! computed once khtr is 3D.
  real                :: heff                 !< Harmonic mean of layer thicknesses           [H ~> m or kg m-2]
  real                :: inv_heff             !< Inverse of the harmonic mean of layer thicknesses
                                              !!  [H-1 ~> m-1 or m2 kg-1]
  real                :: phi_L_avg, phi_R_avg !< Bulk, thickness-weighted tracer averages (left and right column)
                                              !!                                                   [conc m^-3 ]
  real    :: htot                      !< Total column thickness [H ~> m or kg m-2]
  !real    :: heff_tot                  !< Total effective column thickness in the transition layer [m]
  integer :: k, k_bot_min, k_top_max   !< k-indices, min and max for bottom and top, respectively
  integer :: k_bot_max, k_top_min      !< k-indices, max and min for bottom and top, respectively
  integer :: k_bot_diff, k_top_diff    !< different between left and right k-indices for bottom and top, respectively
  integer :: k_top_L, k_bot_L          !< k-indices left native grid
  integer :: k_top_R, k_bot_R          !< k-indices right native grid
  integer :: k_top_zgrid_L, k_bot_zgrid_L !< k-indices left zgrid
  integer :: k_top_zgrid_R, k_bot_zgrid_R !< k-indices right zgrid
  real    :: zeta_top_L, zeta_top_R    !< distance from the top of a layer to the boundary
                                       !! layer depth in the native grid                  [nondim]
  real    :: zeta_bot_L, zeta_bot_R    !< distance from the bottom of a layer to the boundary
                                       !!layer depth in the native grid                   [nondim]
  real    :: zeta_top_zgrid_L, zeta_top_zgrid_R !< distance from the top of a layer to the boundary
                                       !! layer depth in the zgrid                        [nondim]
  real    :: zeta_bot_zgrid_L, zeta_bot_zgrid_R !< distance from the bottom of a layer to the boundary
                                       !!layer depth in the zgrid                         [nondim]
  real    :: h_work_L, h_work_R  !< dummy variables
  real    :: hbl_min             !< minimum BLD (left and right)                          [m]
  real    :: wgt                 !< weight to be used in the linear transition to the interior [nondim]
  real    :: a                   !< coefficient to be used in the linear transition to the interior [nondim]

  F_layer(:) = 0.0
  if (hbl_L == 0. .or. hbl_R == 0.) then
    return
  endif

  ! Calculate vertical indices containing the boundary layer in the zgrid
  call boundary_k_range(boundary, nk, h_L, hbl_L, k_top_L, zeta_top_L, k_bot_L, zeta_bot_L)
  call boundary_k_range(boundary, nk, h_R, hbl_R, k_top_R, zeta_top_R, k_bot_R, zeta_bot_R)
  ! Calculate vertical indices containing the boundary layer in dz_top
  call boundary_k_range(boundary, nk_z, CS%dz_top, hbl_L, k_top_zgrid_L, zeta_top_zgrid_L, k_bot_zgrid_L, zeta_bot_zgrid_L)
  call boundary_k_range(boundary, nk_z, CS%dz_top, hbl_R, k_top_zgrid_R, zeta_top_zgrid_R, k_bot_zgrid_R, zeta_bot_zgrid_R)

  call remapping_core_h(CS%remap_cs, nk, h_L, phi_L, nk_z, CS%dz_top, phi_L_local)
  call remapping_core_h(CS%remap_cs, nk, h_R, phi_R, nk_z, CS%dz_top, phi_R_local)

  if (boundary == SURFACE) then
    k_bot_min = MIN(k_bot_zgrid_L, k_bot_zgrid_R)
    k_bot_max = MAX(k_bot_zgrid_L, k_bot_zgrid_R)
    k_bot_diff = (k_bot_max - k_bot_min)

    ! make sure left and right k indices span same range
    if (k_bot_min .ne. k_bot_zgrid_L) then
      k_bot_zgrid_L = k_bot_min
      zeta_bot_zgrid_L = 1.0
    endif
    if (k_bot_min .ne. k_bot_zgrid_R) then
      k_bot_zgrid_R= k_bot_min
      zeta_bot_zgrid_R = 1.0
    endif

    h_work_L = (CS%dz_top(k_bot_zgrid_L) * zeta_bot_zgrid_L)
    h_work_R = (CS%dz_top(k_bot_zgrid_R) * zeta_bot_zgrid_R)

    ! GMM, the following needs to be modified. We need to calculate ppoly0_E_L and ppoly0_coefs_L here...
    !phi_L_avg = average_value_ppoly( nk_z, phi_L_local, ppoly0_E_L, ppoly0_coefs_L, method, k_bot_L, 0., zeta_bot_L)
    !phi_R_avg = average_value_ppoly( nk_z, phi_R_local, ppoly0_E_R, ppoly0_coefs_R, method, k_bot_R, 0., zeta_bot_R)
    !heff = harmonic_mean(h_work_L, h_work_R)

    ! tracer flux where the minimum BLD intersets layer
    ! GMM, khtr_avg should be computed once khtr is 3D
    if ((CS%linear) .and. (k_bot_diff .gt. 1)) then
      ! apply linear decay at the base of hbl
      do k = k_bot_min-1,1,-1
        !heff = harmonic_mean(h_L(k), h_R(k))
        F_layer(k) = -(CS%dz_top(k) * khtr_u) * (phi_R_local(k) - phi_L_local(k))
      enddo
      htot = 0.0
      do k = k_bot_min+1,k_bot_max, 1
        htot = htot + CS%dz_top(k)
      enddo

      a = -1.0/htot
      htot = 0.0
      do k = k_bot_min,k_bot_max, 1
        !heff = harmonic_mean(h_L(k), h_R(k))
        wgt = (a*(htot + (CS%dz_top(k) * 0.5))) + 1.0
        F_layer(k) = -(CS%dz_top(k) * khtr_u) * (phi_R_local(k) - phi_L_local(k)) * wgt
        htot = htot + CS%dz_top(k)
      enddo
    else
      F_layer(k_bot_min) = -(heff * khtr_u) * (phi_R_avg - phi_L_avg)
      do k = k_bot_min-1,1,-1
        !heff = harmonic_mean(h_L(k), h_R(k))
        F_layer(k) = -(CS%dz_top(k) * khtr_u) * (phi_R_local(k) - phi_L_local(k))
      enddo
    endif
  endif

  if (boundary == BOTTOM) then
    ! TODO: GMM add option to apply linear decay
    k_top_max = MAX(k_top_L, k_top_R)
    ! make sure left and right k indices span same range
    if (k_top_max .ne. k_top_L) then
      k_top_L = k_top_max
      zeta_top_L = 1.0
    endif
    if (k_top_max .ne. k_top_R) then
      k_top_R= k_top_max
      zeta_top_R = 1.0
    endif

    h_work_L = (h_L(k_top_L) * zeta_top_L)
    h_work_R = (h_R(k_top_R) * zeta_top_R)

    phi_L_avg = average_value_ppoly( nk, phi_L, ppoly0_E_L, ppoly0_coefs_L, method, k_top_L, 1.0-zeta_top_L, 1.0)
    phi_R_avg = average_value_ppoly( nk, phi_R, ppoly0_E_R, ppoly0_coefs_R, method, k_top_R, 1.0-zeta_top_R, 1.0)
    heff = harmonic_mean(h_work_L, h_work_R)

    ! tracer flux where the minimum BLD intersets layer
    F_layer(k_top_max) = (-heff * khtr_u) * (phi_R_avg - phi_L_avg)

    do k = k_top_max+1,nk
      heff = harmonic_mean(h_L(k), h_R(k))
      F_layer(k) = -(heff * khtr_u) * (phi_R(k) - phi_L(k))
    enddo
  endif

end subroutine fluxes_layer_method


!> Unit tests for near-boundary horizontal mixing
logical function near_boundary_unit_tests( verbose )
  logical,               intent(in) :: verbose !< If true, output additional information for debugging unit tests

  ! Local variables
  integer, parameter    :: nk = 2               ! Number of layers
  integer, parameter    :: deg = 1              ! Degree of reconstruction (linear here)
  integer, parameter    :: method = 1           ! Method used for integrating polynomials
  real, dimension(nk)   :: phi_L, phi_R         ! Tracer values (left and right column)             [ nondim m^-3 ]
  real, dimension(nk)   :: phi_L_avg, phi_R_avg ! Bulk, thickness-weighted tracer averages (left and right column)
  real, dimension(nk,deg+1) :: phi_pp_L, phi_pp_R   ! Coefficients for the linear pseudo-reconstructions
                                                !                                                   [ nondim m^-3 ]

  real, dimension(nk,2) :: ppoly0_E_L, ppoly0_E_R! Polynomial edge values (left and right)          [concentration]
  real, dimension(nk)   :: h_L, h_R             ! Layer thickness (left and right)                  [m]
  real                  :: khtr_u               ! Horizontal diffusivities at U-point               [m^2 s^-1]
  real                  :: hbl_L, hbl_R         ! Depth of the boundary layer (left and right)      [m]
  real                  :: F_bulk               ! Total diffusive flux across the U point           [nondim s^-1]
  real, dimension(nk)   :: F_layer              ! Diffusive flux within each layer at U-point       [nondim s^-1]
  real                  :: h_u, hblt_u          ! Thickness at the u-point                          [m]
  real                  :: khtr_avg             ! Thickness-weighted diffusivity at the u-point     [m^2 s^-1]
  real                  :: heff                 ! Harmonic mean of layer thicknesses                [m]
  real                  :: inv_heff             ! Inverse of the harmonic mean of layer thicknesses [m^[-1]
  character(len=120)    :: test_name            ! Title of the unit test
  integer               :: k_top                ! Index of cell containing top of boundary
  real                  :: zeta_top             ! Nondimension position
  integer               :: k_bot                ! Index of cell containing bottom of boundary
  real                  :: zeta_bot             ! Nondimension position
  real                  :: area_L,area_R        ! Area of grid cell [m^2]
  area_L = 1.; area_R = 1. ! Set to unity for all unit tests

  near_boundary_unit_tests = .false.

  ! Unit tests for boundary_k_range
  test_name = 'Surface boundary spans the entire top cell'
  h_L = (/5.,5./)
  call boundary_k_range(SURFACE, nk, h_L, 5., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 1, 1., test_name, verbose)

  test_name = 'Surface boundary spans the entire column'
  h_L = (/5.,5./)
  call boundary_k_range(SURFACE, nk, h_L, 10., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 1., test_name, verbose)

  test_name = 'Bottom boundary spans the entire bottom cell'
  h_L = (/5.,5./)
  call boundary_k_range(BOTTOM, nk, h_L, 5., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 2, 1., 2, 0., test_name, verbose)

  test_name = 'Bottom boundary spans the entire column'
  h_L = (/5.,5./)
  call boundary_k_range(BOTTOM, nk, h_L, 10., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 1., 2, 0., test_name, verbose)

  test_name = 'Surface boundary intersects second layer'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 17.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 0.75, test_name, verbose)

  test_name = 'Surface boundary intersects first layer'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 2.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 1, 0.25, test_name, verbose)

  test_name = 'Surface boundary is deeper than column thickness'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 21.0, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 1., test_name, verbose)

  test_name = 'Bottom boundary intersects first layer'
  h_L = (/10.,10./)
  call boundary_k_range(BOTTOM, nk, h_L, 17.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0.75, 2, 0., test_name, verbose)

  test_name = 'Bottom boundary intersects second layer'
  h_L = (/10.,10./)
  call boundary_k_range(BOTTOM, nk, h_L, 2.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 2, 0.25, 2, 0., test_name, verbose)

  ! All cases in this section have hbl which are equal to the column thicknesses
  test_name = 'Equal hbl and same layer thicknesses (gradient from right to left)'
  hbl_L = 10; hbl_R = 10
  h_L = (/5.,5./) ; h_R = (/5.,5./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 1.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  ! Without limiter
  !call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R, &
  !                                  ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/-5.0,-5.0/) )

  !! same as above, but with limiter
  !call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R, &
  !                                  ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer, .true.)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/-1.0,-1.0/) )

  test_name = 'Equal hbl and same layer thicknesses (gradient from left to right)'
  hbl_L = 10.; hbl_R = 10.
  h_L = (/5.,5./) ; h_R = (/5.,5./)
  phi_L = (/1.,1./) ; phi_R = (/0.,0./)
  phi_pp_L(1,1) = 1.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 1.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 0.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 0.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 1.; ppoly0_E_L(1,2) = 1.
  ppoly0_E_L(2,1) = 1.; ppoly0_E_L(2,2) = 1.
  ppoly0_E_R(1,1) = 0.; ppoly0_E_R(1,2) = 0.
  ppoly0_E_R(2,1) = 0.; ppoly0_E_R(2,2) = 0.
  khtr_u = 1.
  !call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
  !                                  ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/5.0,5.0/) )

  test_name = 'Equal hbl and same layer thicknesses (no gradient)'
  hbl_L = 10; hbl_R = 10
  h_L = (/5.,5./) ; h_R = (/5.,5./)
  phi_L = (/1.,1./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 1.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 1.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 1.; ppoly0_E_L(1,2) = 1.
  ppoly0_E_L(2,1) = 1.; ppoly0_E_L(2,2) = 1.
  ppoly0_E_R(1,1) = 1.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  !call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
  !                                  ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/0.0,0.0/) )

  test_name = 'Equal hbl and different layer thicknesses (gradient right to left)'
  hbl_L = 16.; hbl_R = 16.
  h_L = (/10.,6./) ; h_R = (/6.,10./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 1.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  !call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
  !                                  ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/-8.0,-8.0/) )

  test_name = 'Equal hbl and same layer thicknesses (diagonal tracer values)'
  hbl_L = 10.; hbl_R = 10.
  h_L = (/5.,5./) ; h_R = (/5.,5./)
  phi_L = (/1.,0./) ; phi_R = (/0.,1./)
  phi_pp_L(1,1) = 1.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 0.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 1.; ppoly0_E_L(1,2) = 1.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 0.; ppoly0_E_R(1,2) = 0.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  !call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
  !                                  ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/0.0,0.0/) )

  test_name = 'Different hbl and different layer thicknesses (gradient from right to left)'
  hbl_L = 12; hbl_R = 20
  h_L = (/6.,6./) ; h_R = (/10.,10./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 1.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  !call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
  !                                  ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/-7.5,-7.5/) )

  ! Cases where hbl < column thickness (polynomial coefficients specified for pseudo-linear reconstruction)

  test_name = 'hbl < column thickness, hbl same, constant concentration each column'
  hbl_L = 2; hbl_R = 2
  h_L = (/1.,2./) ; h_R = (/1.,2./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 1.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  !call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
  !                                  ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/-1.,-1./) )

  test_name = 'hbl < column thickness, hbl same, linear profile right'
  hbl_L = 2; hbl_R = 2
  h_L = (/1.,2./) ; h_R = (/1.,2./)
  phi_L = (/0.,0./) ; phi_R = (/0.5,2./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 0.; phi_pp_R(1,2) = 1.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 2.
  khtr_u = 1.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 0.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 3.
  !call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
  !                                  ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/-1.,-1./) )

  test_name = 'hbl < column thickness, hbl same, linear profile right, khtr=2'
  hbl_L = 2; hbl_R = 2
  h_L = (/1.,2./) ; h_R = (/1.,2./)
  phi_L = (/0.,0./) ; phi_R = (/0.5,2./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 0.; phi_pp_R(1,2) = 1.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 2.
  khtr_u = 2.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 0.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 3.
  !call fluxes_layer_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, &
  !                                  phi_pp_R, ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_layer)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/-1.,-3./) )

  ! unit tests for layer by layer method
  test_name = 'Different hbl and different column thicknesses (gradient from right to left)'
  hbl_L = 12; hbl_R = 20
  h_L = (/6.,6./) ; h_R = (/10.,10./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 1.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  !call fluxes_layer_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, &
  !                                  phi_pp_R, ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_layer)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/-7.5,-7.5/) )

  test_name = 'Different hbl and different column thicknesses (linear profile right)'

  hbl_L = 15; hbl_R = 6
  h_L = (/10.,10./) ; h_R = (/12.,10./)
  phi_L = (/0.,0./) ; phi_R = (/1.,3./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 0.; phi_pp_R(1,2) = 2.
  phi_pp_R(2,1) = 2.; phi_pp_R(2,2) = 2.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 0.; ppoly0_E_R(1,2) = 2.
  ppoly0_E_R(2,1) = 2.; ppoly0_E_R(2,2) = 4.
  khtr_u = 1.
  !call fluxes_layer_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, &
  !                                  phi_pp_R, ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_layer)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/-3.75,0.0/) )
end function near_boundary_unit_tests

!> Returns true if output of near-boundary unit tests does not match correct computed values
!! and conditionally writes results to stream
logical function test_layer_fluxes(verbose, nk, test_name, F_calc, F_ans)
  logical,                    intent(in) :: verbose   !< If true, write results to stdout
  character(len=80),          intent(in) :: test_name !< Brief description of the unit test
  integer,                    intent(in) :: nk        !< Number of layers
  real, dimension(nk),        intent(in) :: F_calc    !< Fluxes of the unitless tracer from the algorithm [s^-1]
  real, dimension(nk),        intent(in) :: F_ans     !< Fluxes of the unitless tracer calculated by hand [s^-1]
  ! Local variables
  integer :: k
  integer, parameter :: stdunit = stdout

  test_layer_fluxes = .false.
  do k=1,nk
    if ( F_calc(k) /= F_ans(k) ) then
      test_layer_fluxes = .true.
      write(stdunit,*) "MOM_lateral_boundary_diffusion, UNIT TEST FAILED: ", test_name
      write(stdunit,10) k, F_calc(k), F_ans(k)
    elseif (verbose) then
      write(stdunit,10) k, F_calc(k), F_ans(k)
    endif
  enddo

10 format("Layer=",i3," F_calc=",f20.16," F_ans",f20.16)
end function test_layer_fluxes

!> Return true if output of unit tests for boundary_k_range does not match answers
logical function test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, k_top_ans, zeta_top_ans,&
                                       k_bot_ans, zeta_bot_ans, test_name, verbose)
  integer :: k_top               !< Index of cell containing top of boundary
  real    :: zeta_top            !< Nondimension position
  integer :: k_bot               !< Index of cell containing bottom of boundary
  real    :: zeta_bot            !< Nondimension position
  integer :: k_top_ans           !< Index of cell containing top of boundary
  real    :: zeta_top_ans        !< Nondimension position
  integer :: k_bot_ans           !< Index of cell containing bottom of boundary
  real    :: zeta_bot_ans        !< Nondimension position
  character(len=80) :: test_name !< Name of the unit test
  logical :: verbose             !< If true always print output

  integer, parameter :: stdunit = stdout

  test_boundary_k_range = k_top .ne. k_top_ans
  test_boundary_k_range = test_boundary_k_range .or. (zeta_top .ne. zeta_top_ans)
  test_boundary_k_range = test_boundary_k_range .or. (k_bot .ne. k_bot_ans)
  test_boundary_k_range = test_boundary_k_range .or. (zeta_bot .ne. zeta_bot_ans)

  if (test_boundary_k_range) write(stdunit,*) "UNIT TEST FAILED: ", test_name
  if (test_boundary_k_range .or. verbose) then
    write(stdunit,20) "k_top", k_top, "k_top_ans", k_top_ans
    write(stdunit,20) "k_bot", k_bot, "k_bot_ans", k_bot_ans
    write(stdunit,30) "zeta_top", zeta_top, "zeta_top_ans", zeta_top_ans
    write(stdunit,30) "zeta_bot", zeta_bot, "zeta_bot_ans", zeta_bot_ans
  endif

  20 format(A,"=",i3,X,A,"=",i3)
  30 format(A,"=",f20.16,X,A,"=",f20.16)


end function test_boundary_k_range

!> \namespace mom_lateral_boundary_diffusion
!!
!! \section section_LBD The Lateral Boundary Diffusion (LBD) framework
!!
!! The LBD framework accounts for the effects of diabatic mesoscale fluxes
!! within surface and bottom boundary layers. Unlike the equivalent adiabatic
!! fluxes, which is applied along neutral density surfaces, LBD is purely
!! horizontal.
!!
!! The bottom boundary layer fluxes remain to be implemented, although most
!! of the steps needed to do so have already been added and tested.
!!
!! Boundary lateral diffusion can be applied using one of the three methods:
!!
!! * [Method #1: Along layer](@ref section_method2) (default);
!! * [Method #2: Bulk layer](@ref section_method1);
!!
!! A brief summary of these methods is provided below.
!!
!! \subsection section_method1 Along layer approach (Method #1)
!!
!! This is the recommended and more straight forward method where diffusion is
!! applied layer by layer using only information from neighboring cells.
!!
!! Step #1: compute vertical indices containing boundary layer (boundary_k_range).
!! For the TOP boundary layer, these are:
!!
!! k_top, k_bot, zeta_top, zeta_bot
!!
!! Step #2: calculate the diffusive flux at each layer:
!!
!! \f[ F_{k} = -KHTR \times h_{eff}(k) \times (\phi_R(k) - \phi_L(k)),  \f]
!! where h_eff is the [harmonic mean](@ref section_harmonic_mean) of the layer thickness
!! in the left and right columns. This method does not require a limiter since KHTR
!! is already limted based on a diffusive CFL condition prior to the call of this
!! module.
!!
!! Step #3: option to linearly decay the flux from k_bot_min to k_bot_max:
!!
!! If LBD_LINEAR_TRANSITION = True and k_bot_diff > 1, the diffusive flux will decay
!! linearly between the top interface of the layer containing the minimum boundary
!! layer depth (k_bot_min) and the lower interface of the layer containing the
!! maximum layer depth (k_bot_max).
!!
!! \subsection section_method2 Bulk layer approach (Method #2)
!!
!! Apply the lateral boundary diffusive fluxes calculated from a 'bulk model'.This
!! is a lower order representation (Kraus-Turner like approach) which assumes that
!! eddies are acting along well mixed layers (i.e., eddies do not know care about
!! vertical tracer gradients within the boundary layer).
!!
!! Step #1: compute vertical indices containing boundary layer (boundary_k_range).
!! For the TOP boundary layer, these are:
!!
!! k_top, k_bot, zeta_top, zeta_bot
!!
!! Step #2: compute bulk averages (thickness weighted) tracer averages (phi_L and phi_R),
!! then calculate the bulk diffusive flux (F_{bulk}):
!!
!! \f[ F_{bulk} = -KHTR \times h_{eff} \times (\phi_R - \phi_L),  \f]
!! where h_eff is the [harmonic mean](@ref section_harmonic_mean) of the boundary layer depth
!! in the left and right columns (\f[ HBL_L \f] and \f[ HBL_R \f], respectively).
!!
!! Step #3: decompose F_bulk onto individual layers:
!!
!! \f[ F_{layer}(k) = F_{bulk} \times h_{frac}(k) ,  \f]
!!
!! where h_{frac} is
!!
!! \f[ h_{frac}(k) = h_u(k) \times \frac{1}{\sum(h_u)}.  \f]
!!
!! h_u is the [harmonic mean](@ref section_harmonic_mean) of thicknesses at each layer.
!! Special care (layer reconstruction) must be taken at k_min = min(k_botL, k_bot_R).
!!
!! Step #4: option to linearly decay the flux from k_bot_min to k_bot_max:
!!
!! If LBD_LINEAR_TRANSITION = True and k_bot_diff > 1, the diffusive flux will decay
!! linearly between the top interface of the layer containing the minimum boundary
!! layer depth (k_bot_min) and the lower interface of the layer containing the
!! maximum layer depth (k_bot_max).
!!
!! Step #5: limit the tracer flux so that 1) only down-gradient fluxes are applied,
!! and 2) the flux cannot be larger than F_max, which is defined using the tracer
!! gradient:
!!
!! \f[ F_{max} = -0.2 \times [(V_R(k) \times \phi_R(k)) - (V_L(k) \times \phi_L(k))],  \f]
!! where V is the cell volume. Why 0.2?
!!          t=0         t=inf
!!           0           .2
!!         0 1 0       .2.2.2
!!           0           .2
!!
!! \subsection section_harmonic_mean Harmonic Mean
!!
!! The harmonic mean (HM) betwen h1 and h2 is defined as:
!!
!! \f[ HM = \frac{2 \times h1 \times h2}{h1 + h2} \f]
!!
end module MOM_lateral_boundary_diffusion
