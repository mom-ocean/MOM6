!> Routines used to calculate the opacity of the ocean.
module MOM_opacity

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_ptr, post_data
use MOM_diag_mediator, only : query_averaging_enabled, register_diag_field
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_string_functions, only : uppercase
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public set_opacity, opacity_init, opacity_end, opacity_manizza, opacity_morel
public extract_optics_slice, extract_optics_fields, optics_nbands
public absorbRemainingSW, sumSWoverBands

!> This type is used to store information about ocean optical properties
type, public :: optics_type
  integer :: nbands     !< The number of penetrating bands of SW radiation

  real, allocatable :: opacity_band(:,:,:,:) !< SW optical depth per unit thickness [Z-1 ~> m-1]
                        !! The number of radiation bands is most rapidly varying (first) index.

  real, allocatable :: sw_pen_band(:,:,:) !< shortwave radiation [Q R Z T-1 ~> W m-2]
                        !! at the surface in each of the nbands bands that penetrates beyond the surface.
                        !! The most rapidly varying dimension is the band.

  real, allocatable :: min_wavelength_band(:)
      !< The minimum wavelength in each band of penetrating shortwave radiation [nm]
  real, allocatable :: max_wavelength_band(:)
      !< The maximum wavelength in each band of penetrating shortwave radiation [nm]

  real :: PenSW_flux_absorb !< A heat flux that is small enough to be completely absorbed in the next
                        !! sufficiently thick layer [C H T-1 ~> degC m s-1 or degC kg m-2 s-1].
  real :: PenSW_absorb_Invlen !< The inverse of the thickness that is used to absorb the remaining
  !! shortwave heat flux when it drops below PEN_SW_FLUX_ABSORB [H ~> m or kg m-2].

  !! Lookup tables for Ohlmann solar penetration scheme
  !! These would naturally exist as private module variables but that is prohibited in MOM6
  real :: dlog10chl           !< Chl increment within lookup table
  real :: chl_min             !< Lower bound of Chl in lookup table
  real :: log10chl_min        !< Lower bound of Chl in lookup table
  real :: log10chl_max        !< Upper bound of Chl in lookup table
  real, allocatable, dimension(:) :: a1_lut,&       !< Coefficient for band 1
       &                             a2_lut,&       !< Coefficient for band 2
       &                             b1_lut,&       !< Exponential decay scale for band 1
       &                             b2_lut         !< Exponential decay scale for band 2

  integer :: answer_date  !< The vintage of the order of arithmetic and expressions in the optics
                          !! calculations.  Values below 20190101 recover the answers from the
                          !! end of 2018, while higher values use updated and more robust
                          !! forms of the same expressions.

end type optics_type

!> The control structure with parameters for the MOM_opacity module
type, public :: opacity_CS ; private
  logical :: var_pen_sw      !<   If true, use one of the CHL_A schemes (specified by OPACITY_SCHEME) to
                             !! determine the e-folding depth of incoming shortwave radiation.
  integer :: opacity_scheme  !<   An integer indicating which scheme should be used to translate
                             !! water properties into the opacity (i.e., the e-folding depth) and
                             !! (perhaps) the number of bands of penetrating shortwave radiation to use.
  real :: pen_sw_scale       !<   The vertical absorption e-folding depth of the
                             !! penetrating shortwave radiation [Z ~> m].
  real :: pen_sw_scale_2nd   !<   The vertical absorption e-folding depth of the
                             !! (2nd) penetrating shortwave radiation [Z ~> m].
  real :: SW_1ST_EXP_RATIO   !< Ratio for 1st exp decay in Two Exp decay opacity [nondim]
  real :: pen_sw_frac        !<   The fraction of shortwave radiation that is
                             !! penetrating with a constant e-folding approach [nondim]
  real :: blue_frac          !<   The fraction of the penetrating shortwave
                             !! radiation that is in the blue band [nondim].
  real :: opacity_land_value !< The value to use for opacity over land [Z-1 ~> m-1].
                             !! The default is 10 m-1 - a value for muddy water.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                             !! regulate the timing of diagnostic output.
  logical :: warning_issued  !< A flag that is used to avoid repetitive warnings.

  !>@{ Diagnostic IDs
  integer :: id_sw_pen = -1, id_sw_vis_pen = -1
  integer, allocatable :: id_opacity(:)
  !>@}
end type opacity_CS

!>@{ Coded integers to specify the opacity scheme
integer, parameter :: NO_SCHEME = 0, MANIZZA_05 = 1, MOREL_88 = 2, SINGLE_EXP = 3, DOUBLE_EXP = 4,&
     &                OHLMANN_03 = 5
!>@}

character*(10), parameter :: MANIZZA_05_STRING = "MANIZZA_05" !< String to specify the opacity scheme
character*(10), parameter :: MOREL_88_STRING   = "MOREL_88"   !< String to specify the opacity scheme
character*(10), parameter :: OHLMANN_03_STRING = "OHLMANN_03" !< String to specify the opacity scheme
character*(10), parameter :: SINGLE_EXP_STRING = "SINGLE_EXP" !< String to specify the opacity scheme
character*(10), parameter :: DOUBLE_EXP_STRING = "DOUBLE_EXP" !< String to specify the opacity scheme

contains

!> This sets the opacity of sea water based based on one of several different schemes.
subroutine set_opacity(optics, sw_total, sw_vis_dir, sw_vis_dif, sw_nir_dir, sw_nir_dif, &
                       G, GV, US, CS, chl_2d, chl_3d)
  type(optics_type),       intent(inout) :: optics !< An optics structure that has values
                                                   !! set based on the opacities.
  real, dimension(:,:),    pointer       :: sw_total !< Total shortwave flux into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_vis_dir !< Visible, direct shortwave into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_vis_dif !< Visible, diffuse shortwave into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_nir_dir !< Near-IR, direct shortwave into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_nir_dif !< Near-IR, diffuse shortwave into the ocean [Q R Z T-1 ~> W m-2]
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(opacity_CS)                       :: CS     !< The control structure earlier set up by opacity_init.
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: chl_2d !< Vertically uniform chlorophyll-A concentrations [mg m-3]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: chl_3d !< The chlorophyll-A concentrations of each layer [mg m-3]

  ! Local variables
  integer :: i, j, k, n, is, ie, js, je, nz
  real :: inv_sw_pen_scale  ! The inverse of the e-folding scale [Z-1 ~> m-1].
  real :: Inv_nbands        ! The inverse of the number of bands of penetrating
                            ! shortwave radiation [nondim]
  real :: tmp(SZI_(G),SZJ_(G),SZK_(GV)) ! A 3-d temporary array for diagnosing opacity [Z-1 ~> m-1]
  real :: Pen_SW_tot(SZI_(G),SZJ_(G))   ! The penetrating shortwave radiation
                                        ! summed across all bands [Q R Z T-1 ~> W m-2].
  real :: op_diag_len       ! A tiny lengthscale [Z ~> m] used to remap diagnostics of opacity
                            ! from op to 1/op_diag_len * tanh(op * op_diag_len)
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (present(chl_2d) .or. present(chl_3d)) then
    ! The optical properties are based on chlorophyll concentrations.
    call opacity_from_chl(optics, sw_total, sw_vis_dir, sw_vis_dif, sw_nir_dir, sw_nir_dif, &
                          G, GV, US, CS, chl_2d, chl_3d)
  else ! Use sw e-folding scale set by MOM_input
    if (optics%nbands <= 1) then ; Inv_nbands = 1.0
    else ; Inv_nbands = 1.0 / real(optics%nbands) ; endif

    ! Make sure there is no division by 0.
    inv_sw_pen_scale = 1.0 / max(CS%pen_sw_scale, 0.1*GV%Angstrom_Z, &
                                 GV%dZ_subroundoff)
    if ( CS%Opacity_scheme == DOUBLE_EXP ) then
      !$OMP parallel do default(shared)
      do k=1,nz ; do j=js,je ; do i=is,ie
        optics%opacity_band(1,i,j,k) = inv_sw_pen_scale
        optics%opacity_band(2,i,j,k) = 1.0 / max(CS%pen_sw_scale_2nd, &
             0.1*GV%Angstrom_Z, GV%dZ_subroundoff)
      enddo ; enddo ; enddo
      if (.not.associated(sw_total) .or. (CS%pen_SW_scale <= 0.0)) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie ; do n=1,optics%nbands
          optics%sw_pen_band(n,i,j) = 0.0
        enddo ; enddo ; enddo
      else
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          optics%sw_pen_band(1,i,j) = (CS%SW_1st_EXP_RATIO) * sw_total(i,j)
          optics%sw_pen_band(2,i,j) = (1.-CS%SW_1st_EXP_RATIO) * sw_total(i,j)
        enddo ; enddo
      endif
    else
      do k=1,nz ; do j=js,je ; do i=is,ie  ; do n=1,optics%nbands
        optics%opacity_band(n,i,j,k) = inv_sw_pen_scale
      enddo ; enddo ; enddo ; enddo
      if (.not.associated(sw_total) .or. (CS%pen_SW_scale <= 0.0)) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie ; do n=1,optics%nbands
          optics%sw_pen_band(n,i,j) = 0.0
        enddo ; enddo ; enddo
      else
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie ; do n=1,optics%nbands
          optics%sw_pen_band(n,i,j) = CS%pen_SW_frac * Inv_nbands * sw_total(i,j)
        enddo ; enddo ; enddo
      endif
    endif
  endif

  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_sw_pen > 0) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        Pen_SW_tot(i,j) = 0.0
        do n=1,optics%nbands
          Pen_SW_tot(i,j) = Pen_SW_tot(i,j) + optics%sw_pen_band(n,i,j)
        enddo
      enddo ; enddo
      call post_data(CS%id_sw_pen, Pen_SW_tot, CS%diag)
    endif
    if (CS%id_sw_vis_pen > 0) then
      if (CS%opacity_scheme == MANIZZA_05) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          Pen_SW_tot(i,j) = 0.0
          do n=1,min(optics%nbands,2)
            Pen_SW_tot(i,j) = Pen_SW_tot(i,j) + optics%sw_pen_band(n,i,j)
          enddo
        enddo ; enddo
      else
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          Pen_SW_tot(i,j) = 0.0
          do n=1,optics%nbands
            Pen_SW_tot(i,j) = Pen_SW_tot(i,j) + optics%sw_pen_band(n,i,j)
          enddo
        enddo ; enddo
      endif
      call post_data(CS%id_sw_vis_pen, Pen_SW_tot, CS%diag)
    endif
    do n=1,optics%nbands ; if (CS%id_opacity(n) > 0) then
      op_diag_len = 1.0e-10*US%m_to_Z ! A minimal extinction depth to constrain the range of opacity [Z ~> m]
      !$OMP parallel do default(shared)
      do k=1,nz ; do j=js,je ; do i=is,ie
        ! Remap opacity (op) to 1/L * tanh(op * L) where L is one Angstrom.
        ! This gives a nearly identical value when op << 1/L but allows one to
        ! record the values even at reduced precision when opacity is huge (i.e. opaque).
        tmp(i,j,k) = tanh(op_diag_len * optics%opacity_band(n,i,j,k)) / op_diag_len
      enddo ; enddo ; enddo
      call post_data(CS%id_opacity(n), tmp, CS%diag)
    endif ; enddo
  endif

end subroutine set_opacity


!> This sets the "blue" band opacity based on chlorophyll A concentrations
!! The red portion is lumped into the net heating at the surface.
subroutine opacity_from_chl(optics, sw_total, sw_vis_dir, sw_vis_dif, sw_nir_dir, sw_nir_dif, &
                            G, GV, US, CS, chl_2d, chl_3d)
  type(optics_type),       intent(inout) :: optics !< An optics structure that has values
                                                   !! set based on the opacities.
  real, dimension(:,:),    pointer       :: sw_total !< Total shortwave flux into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_vis_dir !< Visible, direct shortwave into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_vis_dif !< Visible, diffuse shortwave into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_nir_dir !< Near-IR, direct shortwave into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_nir_dif !< Near-IR, diffuse shortwave into the ocean [Q R Z T-1 ~> W m-2]
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(opacity_CS)                       :: CS     !< The control structure.
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: chl_2d !< Vertically uniform chlorophyll-A concentrations [mg m-3]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: chl_3d !< A 3-d field of chlorophyll-A concentrations [mg m-3]

  real :: chl_data(SZI_(G),SZJ_(G)) ! The chlorophyll A concentrations in a layer [mg m-3].
  real :: Inv_nbands        ! The inverse of the number of bands of penetrating
                            ! shortwave radiation [nondim]
  real :: Inv_nbands_nir    ! The inverse of the number of bands of penetrating
                            ! near-infrared radiation [nondim]
  real :: SW_pen_tot        ! The sum across the bands of the penetrating
                            ! shortwave radiation [Q R Z T-1 ~> W m-2].
  real :: SW_vis_tot        ! The sum across the visible bands of shortwave
                            ! radiation [Q R Z T-1 ~> W m-2].
  real :: SW_nir_tot        ! The sum across the near infrared bands of shortwave
                            ! radiation [Q R Z T-1 ~> W m-2].
  character(len=128) :: mesg
  integer :: i, j, k, n, is, ie, js, je, nz, nbands
  logical :: multiband_vis_input, multiband_nir_input, total_sw_input

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

!   In this model, the Morel (modified) and Manizza (modified) schemes
! use the "blue" band in the parameterizations to determine the e-folding
! depth of the incoming shortwave attenuation. The red portion is lumped
! into the net heating at the surface.
! Adding Ohlmann scheme. Needs sw_total and chl as inputs. Produces 2 penetrating bands.
! This implementation follows that in CESM-POP using a lookup table in log10(chl) space.
! The table is initialized in subroutine init_ohlmann and the coefficients are recovered
! with routines lookup_ohlmann_swpen and lookup_ohlmann_opacity.
! Note that this form treats the IR solar input implicitly: the sum of partioning
! coefficients < 1.0. The remainder is non-penetrating and is deposited in first layer
! irrespective of thickness. The Ohlmann (2003) paper states that the scheme is not valid
! for vertcal grids with first layer thickness < 2.0 meters.
!
! Ohlmann, J.C. Ocean radiant heating in climate models. J. Climate, 16, 1337-1351, 2003.
!
! Morel, A., Optical modeling of the upper ocean in relation to its biogenous
!   matter content (case-i waters)., J. Geo. Res., {93}, 10,749--10,768, 1988.
!
! Manizza, M., C. L. Quere, A. Watson, and E. T. Buitenhuis, Bio-optical
!   feedbacks among phytoplankton, upper ocean physics and sea-ice in a
!   global model, Geophys. Res. Let., , L05,603, 2005.

  nbands = optics%nbands

  if (nbands <= 1) then ; Inv_nbands = 1.0
  else ; Inv_nbands = 1.0 / real(nbands) ; endif

  if (nbands <= 2) then ; Inv_nbands_nir = 0.0
  else ; Inv_nbands_nir = 1.0 / real(nbands - 2.0) ; endif

  if (.not.(associated(sw_total) .or. (associated(sw_vis_dir) .and. associated(sw_vis_dif) .and. &
                                       associated(sw_nir_dir) .and. associated(sw_nir_dif)) )) then
    if (.not.CS%warning_issued) then
      call MOM_error(WARNING, &
                     "opacity_from_chl called without any shortwave flux arrays allocated.\n"//&
                     "Consider setting PEN_SW_NBANDS = 0 if no shortwave fluxes are being used.")
    endif
    CS%warning_issued = .true.
  endif

  multiband_vis_input = (associated(sw_vis_dir) .and. associated(sw_vis_dif))
  multiband_nir_input = (associated(sw_nir_dir) .and. associated(sw_nir_dif))
  total_sw_input = associated(sw_total)

  chl_data(:,:) = 0.0
  if (present(chl_3d)) then
    do j=js,je ; do i=is,ie ; chl_data(i,j) = chl_3d(i,j,1) ; enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      if ((G%mask2dT(i,j) > 0.0) .and. (chl_3d(i,j,k) < 0.0)) then
        write(mesg,'(" Negative chl_3d of ",(1pe12.4)," found at i,j,k = ", &
                  & 3(1x,i3), " lon/lat = ",(1pe12.4)," E ", (1pe12.4), " N.")') &
                   chl_3d(i,j,k), i, j, k, G%geoLonT(i,j), G%geoLatT(i,j)
        call MOM_error(FATAL, "MOM_opacity opacity_from_chl: "//trim(mesg))
      endif
    enddo ; enddo ; enddo
  elseif (present(chl_2d)) then
    do j=js,je ; do i=is,ie ; chl_data(i,j) = chl_2d(i,j) ; enddo ; enddo
    do j=js,je ; do i=is,ie
      if ((G%mask2dT(i,j) > 0.0) .and. (chl_2d(i,j) < 0.0)) then
        write(mesg,'(" Negative chl_2d of ",(1pe12.4)," at i,j = ", &
                  & 2(i3), "lon/lat = ",(1pe12.4)," E ", (1pe12.4), " N.")') &
                   chl_data(i,j), i, j, G%geoLonT(i,j), G%geoLatT(i,j)
        call MOM_error(FATAL, "MOM_opacity opacity_from_chl: "//trim(mesg))
      endif
    enddo ; enddo
  else
    call MOM_error(FATAL, "Either chl_2d or chl_3d must be present in a call to opacity_form_chl.")
  endif

  select case (CS%opacity_scheme)
    case (MANIZZA_05)
      !$OMP parallel do default(shared) private(SW_vis_tot,SW_nir_tot)
      do j=js,je ; do i=is,ie
        SW_vis_tot = 0.0 ; SW_nir_tot = 0.0
        if (G%mask2dT(i,j) > 0.0) then
          if (multiband_vis_input) then
            SW_vis_tot = sw_vis_dir(i,j) + sw_vis_dif(i,j)
          elseif (total_sw_input) then
            ! Follow Manizza 05 in assuming that 42% of SW is visible.
            SW_vis_tot = 0.42 * sw_total(i,j)
          endif
          if (multiband_nir_input) then
            SW_nir_tot = sw_nir_dir(i,j) + sw_nir_dif(i,j)
          elseif (total_sw_input) then
            SW_nir_tot = sw_total(i,j) - SW_vis_tot
          endif
        endif

        ! Band 1 is Manizza blue.
        optics%sw_pen_band(1,i,j) = CS%blue_frac*sw_vis_tot
        ! Band 2 (if used) is Manizza red.
        if (nbands > 1) &
          optics%sw_pen_band(2,i,j) = (1.0-CS%blue_frac)*sw_vis_tot
        ! All remaining bands are NIR, for lack of something better to do.
        do n=3,nbands
          optics%sw_pen_band(n,i,j) = Inv_nbands_nir * sw_nir_tot
        enddo
      enddo ; enddo
    case (MOREL_88)
      !$OMP parallel do default(shared) private(SW_pen_tot)
      do j=js,je ; do i=is,ie
        SW_pen_tot = 0.0
        if (G%mask2dT(i,j) > 0.0) then
          if (multiband_vis_input) then
            SW_pen_tot = SW_pen_frac_morel(chl_data(i,j)) * (sw_vis_dir(i,j) + sw_vis_dif(i,j))
          elseif (total_sw_input) then
            SW_pen_tot = SW_pen_frac_morel(chl_data(i,j)) * 0.5*sw_total(i,j)
          endif
        endif

        do n=1,nbands
          optics%sw_pen_band(n,i,j) = Inv_nbands*sw_pen_tot
        enddo
     enddo; enddo
    case (OHLMANN_03)
      ! want exactly two penetrating bands. If not, throw an error.
      if ( nbands /= 2 ) then
        call MOM_error(FATAL, "opacity_from_chl: CS%opacity_scheme requires nbands==2.")
      endif
      !$OMP parallel do default(shared) private(SW_vis_tot)
      do j=js,je ; do i=is,ie
        SW_vis_tot = 0.0  ! Ohlmann does not classify as vis/nir. Using vis to add up total
        if (G%mask2dT(i,j) < 0.5) then
          optics%sw_pen_band(1:2,i,j) = 0.  ! Make sure there is a valid value for land points
        else
          if (multiband_vis_input ) then ! If multiband_vis_input is true then so is multiband_nir_input
             SW_vis_tot = sw_vis_dir(i,j) + sw_vis_dif(i,j) + &
                  &       sw_nir_dir(i,j) + sw_nir_dif(i,j)
          elseif (total_sw_input) then
             SW_vis_tot = sw_total(i,j)
          else
            call MOM_error(FATAL, "No shortwave input was provided.")
          endif

          ! Bands 1-2 (Ohlmann factors A with coefficients for Table 1a)
          optics%sw_pen_band(1:2,i,j)  = lookup_ohlmann_swpen(chl_data(i,j),optics)*SW_vis_tot
        endif
      enddo; enddo
    case default
      call MOM_error(FATAL, "opacity_from_chl: CS%opacity_scheme is not valid.")
  end select

  !$OMP parallel do default(shared) firstprivate(chl_data)
  do k=1,nz
    !! FOB
    !!! I don't think this is what we want to do with Ohlmann.
    !!! The surface CHL is used in developing the parameterization.
    !!! Only the surface CHL is used above in setting optics%sw_pen_band for all schemes.
    !!! Seems inconsistent to use depth dependent CHL in opacity calculation.
    !!! Nevertheless, leaving as is for now.
    !! FOB
    if (present(chl_3d)) then
      do j=js,je ; do i=is,ie ; chl_data(i,j) = chl_3d(i,j,k) ; enddo ; enddo
    endif

    select case (CS%opacity_scheme)
      case (MANIZZA_05)
        do j=js,je ; do i=is,ie
          if (G%mask2dT(i,j) <= 0.5) then
            do n=1,optics%nbands
              optics%opacity_band(n,i,j,k) = CS%opacity_land_value
            enddo
          else
            ! Band 1 is Manizza blue.
            optics%opacity_band(1,i,j,k) = (0.0232 + 0.074*chl_data(i,j)**0.674) * US%Z_to_m
            if (nbands >= 2) &  !  Band 2 is Manizza red.
              optics%opacity_band(2,i,j,k) = (0.225 + 0.037*chl_data(i,j)**0.629) * US%Z_to_m
            ! All remaining bands are NIR, for lack of something better to do.
            do n=3,nbands ; optics%opacity_band(n,i,j,k) = 2.86*US%Z_to_m ; enddo
          endif
        enddo ; enddo
      case (MOREL_88)
        do j=js,je ; do i=is,ie
          optics%opacity_band(1,i,j,k) = CS%opacity_land_value
          if (G%mask2dT(i,j) > 0.0) &
            optics%opacity_band(1,i,j,k) = US%Z_to_m * opacity_morel(chl_data(i,j))

          do n=2,optics%nbands
            optics%opacity_band(n,i,j,k) = optics%opacity_band(1,i,j,k)
          enddo
        enddo; enddo
      case (OHLMANN_03)
        !! not testing for 2 bands since we did it above
        do j=js,je ; do i=is,ie
          if (G%mask2dT(i,j) <= 0.5) then
            optics%opacity_band(1:2,i,j,k) = CS%opacity_land_value
          else
            ! Bands 1-2 (Ohlmann factors B with coefficients for Table 1a
            optics%opacity_band(1:2,i,j,k) = lookup_ohlmann_opacity(chl_data(i,j),optics) * US%Z_to_m
          endif
        enddo; enddo
      case default
        call MOM_error(FATAL, "opacity_from_chl: CS%opacity_scheme is not valid.")
    end select
  enddo

end subroutine opacity_from_chl

!> This sets the blue-wavelength opacity according to the scheme proposed by
!! Morel and Antoine (1994).
function opacity_morel(chl_data)
  real, intent(in)  :: chl_data !< The chlorophyll-A concentration in [mg m-3]
  real :: opacity_morel !< The returned opacity [m-1]

  !   The following are coefficients for the optical model taken from Morel and
  ! Antoine (1994). These coefficients represent a non uniform distribution of
  ! chlorophyll-a through the water column.  Other approaches may be more
  ! appropriate when using an interactive ecosystem model that predicts
  ! three-dimensional chl-a values.
  real, dimension(6), parameter :: &
    Z2_coef = (/7.925, -6.644, 3.662, -1.815, -0.218,  0.502/) ! Extinction length coefficients [m]
  real :: Chl, Chl2 ! The log10 of chl_data (in mg m-3), and Chl^2 [nondim]

  Chl = log10(min(max(chl_data,0.02),60.0)) ; Chl2 = Chl*Chl
  opacity_morel = 1.0 / ( (Z2_coef(1) + Z2_coef(2)*Chl) + Chl2 * &
      ((Z2_coef(3) + Chl*Z2_coef(4)) + Chl2*(Z2_coef(5) + Chl*Z2_coef(6))) )
end function

!> This sets the penetrating shortwave fraction according to the scheme proposed by
!! Morel and Antoine (1994).
function SW_pen_frac_morel(chl_data)
  real, intent(in)  :: chl_data !< The chlorophyll-A concentration [mg m-3]
  real :: SW_pen_frac_morel     !< The returned penetrating shortwave fraction [nondim]

  !   The following are coefficients for the optical model taken from Morel and
  ! Antoine (1994). These coefficients represent a non uniform distribution of
  ! chlorophyll-a through the water column.  Other approaches may be more
  ! appropriate when using an interactive ecosystem model that predicts
  ! three-dimensional chl-a values.
  real :: Chl, Chl2         ! The log10 of chl_data in mg m-3, and Chl^2 [nondim]
  real, dimension(6), parameter :: &
    V1_coef = (/0.321,  0.008, 0.132,  0.038, -0.017, -0.007/) ! Penetrating fraction coefficients [nondim]

  Chl = log10(min(max(chl_data,0.02),60.0)) ; Chl2 = Chl*Chl
  SW_pen_frac_morel = 1.0 - ( (V1_coef(1) + V1_coef(2)*Chl) + Chl2 * &
       ((V1_coef(3) + Chl*V1_coef(4)) + Chl2*(V1_coef(5) + Chl*V1_coef(6))) )
end function SW_pen_frac_morel

!>   This sets the blue-wavelength opacity according to the scheme proposed by
!! Manizza, M. et al, 2005.
function opacity_manizza(chl_data)
  real, intent(in)  :: chl_data !< The chlorophyll-A concentration [mg m-3]
  real :: opacity_manizza !< The returned opacity [m-1]
!   This sets the blue-wavelength opacity according to the scheme proposed by Manizza, M. et al, 2005.

  opacity_manizza = 0.0232 + 0.074*chl_data**0.674
end function

!> This subroutine returns a 2-d slice at constant j of fields from an optics_type, with the potential
!! for rescaling these fields.
subroutine extract_optics_slice(optics, j, G, GV, opacity, opacity_scale, penSW_top, penSW_scale, SpV_avg)
  type(optics_type),       intent(in)  :: optics !< An optics structure that has values of opacities
                                                 !! and shortwave fluxes.
  integer,                 intent(in)  :: j      !< j-index to extract
  type(ocean_grid_type),   intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV     !< The ocean's vertical grid structure.
  real, dimension(max(optics%nbands,1),SZI_(G),SZK_(GV)), &
                 optional, intent(out) :: opacity   !< The opacity in each band, i-point, and layer [Z-1 ~> m-1],
                                                    !! but with units that can be altered by opacity_scale
                                                    !! and the presence of SpV_avg to change this to other
                                                    !! units like [H-1 ~> m-1 or m2 kg-1]
  real,          optional, intent(in)  :: opacity_scale !< A factor by which to rescale the opacity [nondim] or
                                                    !! [Z H-1 ~> 1 or m3 kg-1]
  real, dimension(max(optics%nbands,1),SZI_(G)), &
                 optional, intent(out) :: penSW_top !< The shortwave radiation [Q R Z T-1 ~> W m-2]
                                                    !! at the surface in each of the nbands bands
                                                    !! that penetrates beyond the surface skin layer.
  real,          optional, intent(in)  :: penSW_scale !< A factor by which to rescale the shortwave flux [nondim]
                                                    !! or other units.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)  :: SpV_avg   !< The layer-averaged specific volume [R-1 ~> m3 kg-1]
                                                    !! that is used along with opacity_scale in non-Boussinesq
                                                    !! cases to change the opacity from distance based units to
                                                    !! mass-based units

  ! Local variables
  real :: scale_opacity ! A rescaling factor for opacity [nondim], or the same units as opacity_scale.
  real :: scale_penSW   ! A rescaling factor for the penetrating shortwave radiation [nondim] or the
                        ! same units as penSW_scale
  integer :: i, is, ie, k, nz, n
  is = G%isc ; ie = G%iec ; nz = GV%ke

  scale_opacity = 1.0 ; if (present(opacity_scale)) scale_opacity = opacity_scale
  scale_penSW = 1.0 ; if (present(penSW_scale)) scale_penSW = penSW_scale

  if (present(opacity)) then
    if (present(SpV_avg)) then
      do k=1,nz ; do i=is,ie ; do n=1,optics%nbands
        opacity(n,i,k) = (scale_opacity * SpV_avg(i,j,k)) * optics%opacity_band(n,i,j,k)
      enddo ; enddo ; enddo
    else
      do k=1,nz ; do i=is,ie ; do n=1,optics%nbands
        opacity(n,i,k) = scale_opacity * optics%opacity_band(n,i,j,k)
      enddo ; enddo ; enddo
    endif
  endif

  if (present(penSW_top)) then ; do i=is,ie ; do n=1,optics%nbands
    penSW_top(n,i) = scale_penSW * optics%sw_pen_band(n,i,j)
  enddo ; enddo ; endif

end subroutine extract_optics_slice

!> Set arguments to fields from the optics type.
subroutine extract_optics_fields(optics, nbands)
  type(optics_type),       intent(in)  :: optics !< An optics structure that has values of opacities
                                                 !! and shortwave fluxes.
  integer, optional,       intent(out) :: nbands !< The number of penetrating bands of SW radiation

  if (present(nbands)) nbands = optics%nbands

end subroutine extract_optics_fields

!> Return the number of bands of penetrating shortwave radiation.
function optics_nbands(optics)
  type(optics_type),           pointer :: optics !< An optics structure that has values of opacities
                                                 !! and shortwave fluxes.
  integer :: optics_nbands !< The number of penetrating bands of SW radiation

  if (associated(optics)) then
    optics_nbands = optics%nbands
  else
    optics_nbands = 0
  endif
end function optics_nbands

!> Apply shortwave heating below the boundary layer (when running with the bulk mixed layer inherited
!! from GOLD) or throughout the water column.
!!
!! In addition, it causes all of the remaining SW radiation to be absorbed, provided that the total
!! water column thickness is greater than H_limit_fluxes.
!! For thinner water columns, the heating is scaled down proportionately, the assumption being that the
!! remaining heating (which is left in Pen_SW) should go into an (absent for now) ocean bottom sediment layer.
subroutine absorbRemainingSW(G, GV, US, h, opacity_band, nsw, optics, j, dt, H_limit_fluxes, &
                             adjustAbsorptionProfile, absorbAllSW, T, Pen_SW_bnd, &
                             eps, ksort, htot, Ttot, TKE, dSV_dT)

  type(ocean_grid_type),             intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),           intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),             intent(in)    :: US   !< A dimensional unit scaling type
  integer,                           intent(in)    :: nsw  !< Number of bands of penetrating
                                                           !! shortwave radiation.
  real, dimension(SZI_(G),SZK_(GV)), intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(max(1,nsw),SZI_(G),SZK_(GV)), intent(in) :: opacity_band !< Opacity in each band of penetrating
                                                           !! shortwave radiation [H-1 ~> m-1 or m2 kg-1].
                                                           !! The indices are band, i, k.
  type(optics_type),                 intent(in)    :: optics !< An optics structure that has values of
                                                           !! opacities and shortwave fluxes.
  integer,                           intent(in)    :: j    !< j-index to work on.
  real,                              intent(in)    :: dt   !< Time step [T ~> s].
  real,                              intent(in)    :: H_limit_fluxes !< If the total ocean depth is
                                                           !! less than this, they are scaled away
                                                           !! to avoid numerical instabilities
                                                           !! [H ~> m or kg m-2]. This would
                                                           !! not be necessary if a finite heat
                                                           !! capacity mud-layer were added.
  logical,                          intent(in)    :: adjustAbsorptionProfile !< If true, apply
                                                           !! heating above the layers in which it
                                                           !! should have occurred to get the
                                                           !! correct mean depth (and potential
                                                           !! energy change) of the shortwave that
                                                           !! should be absorbed by each layer.
  logical,                          intent(in)    :: absorbAllSW !< If true, apply heating above the
                                                           !! layers in which it should have occurred
                                                           !! to get the correct mean depth (and
                                                           !! potential energy change) of the
                                                           !! shortwave that should be absorbed by
                                                           !! each layer.
  real, dimension(SZI_(G),SZK_(GV)), intent(inout) :: T    !< Layer potential/conservative
                                                           !! temperatures [C ~> degC]
  real, dimension(max(1,nsw),SZI_(G)), intent(inout) :: Pen_SW_bnd !< Penetrating shortwave heating in
                                                           !! each band that hits the bottom and will
                                                           !! will be redistributed through the water
                                                           !! column [C H ~> degC m or degC kg m-2],
                                                           !! size nsw x SZI_(G).
  real, dimension(SZI_(G),SZK_(GV)), optional, intent(in) :: eps !< Small thickness that must remain in
                                                           !! each layer, and which will not be
                                                           !! subject to heating [H ~> m or kg m-2]
  integer, dimension(SZI_(G),SZK_(GV)), optional, intent(in) :: ksort !< Density-sorted k-indices.
  real, dimension(SZI_(G)), optional, intent(in)    :: htot !< Total mixed layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G)), optional, intent(inout) :: Ttot !< Depth integrated mixed layer
                                                           !! temperature [C H ~> degC m or degC kg m-2]
  real, dimension(SZI_(G),SZK_(GV)), optional, intent(in) :: dSV_dT !< The partial derivative of specific volume
                                                           !! with temperature [R-1 C-1 ~> m3 kg-1 degC-1]
  real, dimension(SZI_(G),SZK_(GV)), optional, intent(inout) :: TKE !< The TKE sink from mixing the heating
                                                           !! throughout a layer [R Z3 T-2 ~> J m-2].

  ! Local variables
  real, dimension(SZI_(G),SZK_(GV)) :: &
    T_chg_above    ! A temperature change that will be applied to all the thick
                   ! layers above a given layer [C ~> degC].  This is only nonzero if
                   ! adjustAbsorptionProfile is true, in which case the net
                   ! change in the temperature of a layer is the sum of the
                   ! direct heating of that layer plus T_chg_above from all of
                   ! the layers below, plus any contribution from absorbing
                   ! radiation that hits the bottom.
  real, dimension(SZI_(G)) :: &
    h_heat, &      ! The thickness of the water column that will be heated by
                   ! any remaining shortwave radiation [H ~> m or kg m-2].
    T_chg, &       ! The temperature change of thick layers due to the remaining
                   ! shortwave radiation and contributions from T_chg_above [C ~> degC].
    Pen_SW_rem     ! The sum across all wavelength bands of the penetrating shortwave
                   ! heating that hits the bottom and will be redistributed through
                   ! the water column [C H ~> degC m or degC kg m-2]
  real :: SW_trans          ! fraction of shortwave radiation that is not
                            ! absorbed in a layer [nondim]
  real :: unabsorbed        ! fraction of the shortwave radiation that
                            ! is not absorbed because the layers are too thin [nondim]
  real :: Ih_limit          ! inverse of the total depth at which the
                            ! surface fluxes start to be limited [H-1 ~> m-1 or m2 kg-1]
  real :: h_min_heat        ! minimum thickness layer that should get heated [H ~> m or kg m-2]
  real :: opt_depth         ! optical depth of a layer [nondim]
  real :: exp_OD            ! exp(-opt_depth) [nondim]
  real :: heat_bnd          ! heating due to absorption in the current
                            ! layer by the current band, including any piece that
                            ! is moved upward [C H ~> degC m or degC kg m-2]
  real :: SWa               ! fraction of the absorbed shortwave that is
                            ! moved to layers above with adjustAbsorptionProfile [nondim]
  real :: coSWa_frac        ! The fraction of SWa that is actually moved upward [nondim]
  real :: min_SW_heat       ! A minimum remaining shortwave heating within a timestep that will be simply
                            ! absorbed in the next layer for computational efficiency, instead of
                            ! continuing to penetrate [C H ~> degC m or degC kg m-2].
  real :: I_Habs            ! The inverse of the absorption length for a minimal flux [H-1 ~> m-1 or m2 kg-1]
  real :: epsilon           ! A small thickness that must remain in each
                            ! layer, and which will not be subject to heating [H ~> m or kg m-2]
  real :: g_Hconv2          ! A conversion factor for use in the TKE calculation
                            ! in units of [Z3 R2 T-2 H-2 ~> kg2 m-5 s-2 or m s-2].
  logical :: SW_Remains     ! If true, some column has shortwave radiation that
                            ! was not entirely absorbed.
  logical :: TKE_calc       ! If true, calculate the implications to the
                            ! TKE budget of the shortwave heating.
  real :: C1_6, C1_60       ! Rational fractions [nondim]
  integer :: is, ie, nz, i, k, ks, n

  if (nsw < 1) return

  SW_Remains = .false.
  min_SW_heat = optics%PenSW_flux_absorb * dt
  I_Habs = optics%PenSW_absorb_Invlen

  h_min_heat = 2.0*GV%Angstrom_H + GV%H_subroundoff
  is = G%isc ; ie = G%iec ; nz = GV%ke
  C1_6 = 1.0 / 6.0 ; C1_60 = 1.0 / 60.0

  TKE_calc = (present(TKE) .and. present(dSV_dT))

  if (optics%answer_date < 20190101) then
    g_Hconv2 = (US%L_to_Z**2*GV%g_Earth * GV%H_to_RZ) * GV%H_to_RZ
  else
    g_Hconv2 = US%L_to_Z**2*GV%g_Earth * GV%H_to_RZ**2
  endif

  h_heat(:) = 0.0
  if (present(htot)) then ; do i=is,ie ; h_heat(i) = htot(i) ; enddo ; endif

  ! Apply penetrating SW radiation to remaining parts of layers.
  ! Excessively thin layers are not heated to avoid runaway temps.
  do ks=1,nz ; do i=is,ie
    k = ks
    if (present(ksort)) then
      if (ksort(i,ks) <= 0) cycle
      k = ksort(i,ks)
    endif
    epsilon = 0.0 ; if (present(eps)) epsilon = eps(i,k)

    T_chg_above(i,k) = 0.0

    if (h(i,k) > 1.5*epsilon) then
      do n=1,nsw ; if (Pen_SW_bnd(n,i) > 0.0) then
        ! SW_trans is the SW that is transmitted THROUGH the layer
        opt_depth = h(i,k) * opacity_band(n,i,k)
        exp_OD = exp(-opt_depth)
        SW_trans = exp_OD

        ! Heating at a very small rate can be absorbed by a sufficiently thick layer or several
        ! thin layers without further penetration.
        if (optics%answer_date < 20190101) then
          if (nsw*Pen_SW_bnd(n,i)*SW_trans < min_SW_heat*min(1.0, I_Habs*h(i,k)) ) SW_trans = 0.0
        elseif ((nsw*Pen_SW_bnd(n,i)*SW_trans < min_SW_heat) .and. (h(i,k) > h_min_heat)) then
          if (nsw*Pen_SW_bnd(n,i) <= min_SW_heat * (I_Habs*(h(i,k) - h_min_heat))) then
            SW_trans = 0.0
          else
            SW_trans = min(SW_trans, &
                           1.0 - (min_SW_heat*(I_Habs*(h(i,k) - h_min_heat))) / (nsw*Pen_SW_bnd(n,i)))
          endif
        endif

        Heat_bnd = Pen_SW_bnd(n,i) * (1.0 - SW_trans)
        if (adjustAbsorptionProfile .and. (h_heat(i) > 0.0)) then
          !   In this case, a fraction of the heating is applied to the
          ! overlying water so that the mean pressure at which the shortwave
          ! heating occurs is exactly what it would have been with a careful
          ! pressure-weighted averaging of the exponential heating profile,
          ! hence there should be no TKE budget requirements due to this
          ! layer.  Very clever, but this is also limited so that the
          ! water above is not heated at a faster rate than the layer
          ! actually being heated, i.e., SWA <= h_heat / (h_heat + h(i,k))
          ! and takes the energetics of the rest of the heating into account.
          ! (-RWH, ~7 years later.)
          if (opt_depth > 1e-5) then
            SWa = ((opt_depth + (opt_depth + 2.0)*exp_OD) - 2.0) / &
              ((opt_depth + opacity_band(n,i,k) * h_heat(i)) * &
               (1.0 - exp_OD))
          else
            ! Use Taylor series expansion of the expression above for a
            ! more accurate form with very small layer optical depths.
            SWa = h(i,k) * (opt_depth * (1.0 - opt_depth)) / &
              ((h_heat(i) + h(i,k)) * (6.0 - 3.0*opt_depth))
          endif
          coSWa_frac = 0.0
          if (SWa*(h_heat(i) + h(i,k)) > h_heat(i)) then
            coSWa_frac = (SWa*(h_heat(i) + h(i,k)) - h_heat(i) ) / &
                         (SWa*(h_heat(i) + h(i,k)))
            SWa = h_heat(i) / (h_heat(i) + h(i,k))
          endif

          T_chg_above(i,k) = T_chg_above(i,k) + (SWa * Heat_bnd) / h_heat(i)
          T(i,k) = T(i,k) + ((1.0 - SWa) * Heat_bnd) / h(i,k)
        else
          coSWa_frac = 1.0
          T(i,k) = T(i,k) + Pen_SW_bnd(n,i) * (1.0 - SW_trans) / h(i,k)
        endif

        if (TKE_calc) then
          if (opt_depth > 1e-2) then
            TKE(i,k) = TKE(i,k) - coSWa_frac*Heat_bnd*dSV_dT(i,k)* &
               (0.5*h(i,k)*g_Hconv2) * &
               (opt_depth*(1.0+exp_OD) - 2.0*(1.0-exp_OD)) / (opt_depth*(1.0-exp_OD))
          else
            ! Use Taylor series-derived approximation to the above expression
            ! that is well behaved and more accurate when opt_depth is small.
            TKE(i,k) = TKE(i,k) - coSWa_frac*Heat_bnd*dSV_dT(i,k)* &
               (0.5*h(i,k)*g_Hconv2) * &
               (C1_6*opt_depth * (1.0 - C1_60*opt_depth**2))
          endif
        endif

        Pen_SW_bnd(n,i) = Pen_SW_bnd(n,i) * SW_trans
      endif ; enddo
    endif

    ! Add to the accumulated thickness above that could be heated.
    ! Only layers greater than h_min_heat thick should get heated.
    if (h(i,k) >= 2.0*h_min_heat) then
      h_heat(i) = h_heat(i) + h(i,k)
    elseif (h(i,k) > h_min_heat) then
      h_heat(i) = h_heat(i) + (2.0*h(i,k) - 2.0*h_min_heat)
    endif
  enddo ; enddo ! i & k loops

! if (.not.absorbAllSW .and. .not.adjustAbsorptionProfile) return

  ! Unless modified, there is no temperature change due to fluxes from the bottom.
  do i=is,ie ; T_chg(i) = 0.0 ; enddo

  if (absorbAllSW) then
    ! If there is still shortwave radiation at this point, it could go into
    ! the bottom (with a bottom mud model), or it could be redistributed back
    ! through the water column.
    do i=is,ie
      Pen_SW_rem(i) = Pen_SW_bnd(1,i)
      do n=2,nsw ; Pen_SW_rem(i) = Pen_SW_rem(i) + Pen_SW_bnd(n,i) ; enddo
    enddo
    do i=is,ie ; if (Pen_SW_rem(i) > 0.0) SW_Remains = .true. ; enddo

    Ih_limit = 1.0 / H_limit_fluxes
    do i=is,ie ; if ((Pen_SW_rem(i) > 0.0) .and. (h_heat(i) > 0.0)) then
      if (h_heat(i)*Ih_limit >= 1.0) then
        T_chg(i) = Pen_SW_rem(i) / h_heat(i) ; unabsorbed = 0.0
      else
        T_chg(i) = Pen_SW_rem(i) * Ih_limit
        unabsorbed = 1.0 - h_heat(i)*Ih_limit
      endif
      do n=1,nsw ; Pen_SW_bnd(n,i) = unabsorbed * Pen_SW_bnd(n,i) ; enddo
    endif ; enddo
  endif ! absorbAllSW

  if (absorbAllSW .or. adjustAbsorptionProfile) then
    do ks=nz,1,-1 ; do i=is,ie
      k = ks
      if (present(ksort)) then
        if (ksort(i,ks) <= 0) cycle
        k = ksort(i,ks)
      endif

      if (T_chg(i) > 0.0) then
        ! Only layers greater than h_min_heat thick should get heated.
        if (h(i,k) >= 2.0*h_min_heat) then ; T(i,k) = T(i,k) + T_chg(i)
        elseif (h(i,k) > h_min_heat) then
          T(i,k) = T(i,k) + T_chg(i) * (2.0 - 2.0*h_min_heat/h(i,k))
        endif
      endif
      ! Increase the heating for layers above.
      T_chg(i) = T_chg(i) + T_chg_above(i,k)
    enddo ; enddo
    if (present(htot) .and. present(Ttot)) then
      do i=is,ie ; Ttot(i) = Ttot(i) + T_chg(i) * htot(i) ; enddo
    endif
  endif ! absorbAllSW .or. adjustAbsorptionProfile

end subroutine absorbRemainingSW


!> This subroutine calculates the total shortwave heat flux integrated over
!! bands as a function of depth.  This routine is only called for computing
!! buoyancy fluxes for use in KPP. This routine does not update the state.
subroutine sumSWoverBands(G, GV, US, h, dz, nsw, optics, j, dt, &
                          H_limit_fluxes, absorbAllSW, iPen_SW_bnd, netPen)
  type(ocean_grid_type),    intent(in)    :: G   !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: h   !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: dz  !< Layer vertical extent [Z ~> m].
  integer,                  intent(in)    :: nsw !< The number of bands of penetrating shortwave
                                                 !! radiation, perhaps from optics_nbands(optics),
  type(optics_type),        intent(in)    :: optics !< An optics structure that has values
                                                   !! set based on the opacities.
  integer,                  intent(in)    :: j   !< j-index to work on.
  real,                     intent(in)    :: dt  !< Time step [T ~> s].
  real,                     intent(in)    :: H_limit_fluxes !< the total depth at which the
                                                 !! surface fluxes start to be limited to avoid
                                                 !! excessive heating of a thin ocean [H ~> m or kg m-2]
  logical,                  intent(in)    :: absorbAllSW !< If true, ensure that all shortwave
                                                 !! radiation is absorbed in the ocean water column.
  real, dimension(max(nsw,1),SZI_(G)), intent(in) :: iPen_SW_bnd !< The incident penetrating shortwave
                                                 !! in each band at the sea surface; size nsw x SZI_(G)
                                                 !! [C H ~> degC m or degC kg m-2].
  real, dimension(SZI_(G),SZK_(GV)+1), &
                             intent(inout) :: netPen !< Net penetrating shortwave heat flux at each
                                                 !! interface, summed across all bands
                                                 !! [C H ~> degC m or degC kg m-2].
  ! Local variables
  real :: h_heat(SZI_(G))     ! thickness of the water column that receives
                              ! remaining shortwave radiation [H ~> m or kg m-2].
  real :: Pen_SW_rem(SZI_(G)) ! sum across all wavelength bands of the
                              ! penetrating shortwave heating that hits the bottom
                              ! and will be redistributed through the water column
                              ! [C H ~> degC m or degC kg m-2]

  real, dimension(max(nsw,1),SZI_(G)) :: Pen_SW_bnd ! The remaining penetrating shortwave radiation
                          ! in each band, initially iPen_SW_bnd [C H ~> degC m or degC kg m-2]
  real :: SW_trans        ! fraction of shortwave radiation not
                          ! absorbed in a layer [nondim]
  real :: unabsorbed      ! fraction of the shortwave radiation
                          ! not absorbed because the layers are too thin [nondim].
  real :: Ih_limit        ! inverse of the total depth at which the
                          ! surface fluxes start to be limited [H-1 ~> m-1 or m2 kg-1]
  real :: min_SW_heat     ! A minimum remaining shortwave heating within a timestep that will be simply
                          ! absorbed in the next layer for computational efficiency, instead of
                          ! continuing to penetrate [C H ~> degC m or degC kg m-2].
  real :: I_Habs            ! The inverse of the absorption length for a minimal flux [H-1 ~> m-1 or m2 kg-1]
  real :: h_min_heat      ! minimum thickness layer that should get heated [H ~> m or kg m-2]
  real :: opt_depth       ! optical depth of a layer [nondim]
  real :: exp_OD          ! exp(-opt_depth) [nondim]
  logical :: SW_Remains   ! If true, some column has shortwave radiation that
                          ! was not entirely absorbed.

  integer :: is, ie, nz, i, k, n
  SW_Remains = .false.

  I_Habs = 1e3*GV%H_to_m ! optics%PenSW_absorb_Invlen

  h_min_heat = 2.0*GV%Angstrom_H + GV%H_subroundoff
  is = G%isc ; ie = G%iec ; nz = GV%ke

  if (nsw < 1) then
    netPen(:,:) = 0.0
    return
  endif

  pen_SW_bnd(:,:) = iPen_SW_bnd(:,:)
  do i=is,ie ; h_heat(i) = 0.0 ; enddo
  do i=is,ie
    netPen(i,1) = 0.
    do n=1,max(nsw,1)
      netPen(i,1) = netPen(i,1) + pen_SW_bnd(n,i)   ! Surface interface
    enddo
  enddo

  ! Apply penetrating SW radiation to remaining parts of layers.
  ! Excessively thin layers are not heated to avoid runaway temps.
  min_SW_heat = optics%PenSW_flux_absorb*dt ! Default of 2.5e-11*US%T_to_s*GV%m_to_H
  do k=1,nz

    do i=is,ie
      netPen(i,k+1) = 0.

      if (h(i,k) > 0.0) then
        do n=1,nsw ; if (Pen_SW_bnd(n,i) > 0.0) then
          ! SW_trans is the SW that is transmitted THROUGH the layer
          opt_depth = dz(i,k) * optics%opacity_band(n,i,j,k)
          exp_OD = exp(-opt_depth)
          SW_trans = exp_OD

          ! Heating at a very small rate can be absorbed by a sufficiently thick layer or several
          ! thin layers without further penetration.
          if (optics%answer_date < 20190101) then
            if (nsw*Pen_SW_bnd(n,i)*SW_trans < min_SW_heat*min(1.0, I_Habs*h(i,k)) ) SW_trans = 0.0
          elseif ((nsw*Pen_SW_bnd(n,i)*SW_trans < min_SW_heat) .and. (h(i,k) > h_min_heat)) then
            if (nsw*Pen_SW_bnd(n,i) <= min_SW_heat * (I_Habs*(h(i,k) - h_min_heat))) then
              SW_trans = 0.0
            else
              SW_trans = min(SW_trans, &
                             1.0 - (min_SW_heat*(I_Habs*(h(i,k) - h_min_heat))) / (nsw*Pen_SW_bnd(n,i)))
            endif
          endif

          Pen_SW_bnd(n,i) = Pen_SW_bnd(n,i) * SW_trans
          netPen(i,k+1)   = netPen(i,k+1) + Pen_SW_bnd(n,i)
        endif ; enddo
      endif ! h(i,k) > 0.0

      ! Add to the accumulated thickness above that could be heated.
      ! Only layers greater than h_min_heat thick should get heated.
      if (h(i,k) >= 2.0*h_min_heat) then
        h_heat(i) = h_heat(i) + h(i,k)
      elseif (h(i,k) > h_min_heat) then
        h_heat(i) = h_heat(i) + (2.0*h(i,k) - 2.0*h_min_heat)
      endif
    enddo ! i loop
  enddo ! k loop

  if (absorbAllSW) then

    ! If there is still shortwave radiation at this point, it could go into
    ! the bottom (with a bottom mud model), or it could be redistributed back
    ! through the water column.
    do i=is,ie
      Pen_SW_rem(i) = Pen_SW_bnd(1,i)
      do n=2,nsw ; Pen_SW_rem(i) = Pen_SW_rem(i) + Pen_SW_bnd(n,i) ; enddo
    enddo
    do i=is,ie ; if (Pen_SW_rem(i) > 0.0) SW_Remains = .true. ; enddo

    Ih_limit = 1.0 / H_limit_fluxes
    do i=is,ie ; if ((Pen_SW_rem(i) > 0.0) .and. (h_heat(i) > 0.0)) then
      if (h_heat(i)*Ih_limit < 1.0) then
        unabsorbed = 1.0 - h_heat(i)*Ih_limit
      else
        unabsorbed = 0.0
      endif
      do n=1,nsw ; Pen_SW_bnd(n,i) = unabsorbed * Pen_SW_bnd(n,i) ; enddo
    endif ; enddo

  endif ! absorbAllSW

end subroutine sumSWoverBands



!> This routine initializes the opacity module, including an optics_type.
subroutine opacity_init(Time, G, GV, US, param_file, diag, CS, optics)
  type(time_type), target, intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< model vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(opacity_CS) :: CS                         !< Opacity control structure
  type(optics_type) :: optics                    !< An optics structure that has parameters
                                                 !! set and arrays allocated here.
  ! Local variables
  character(len=200) :: tmpstr
  character(len=40)  :: mdl = "MOM_opacity"
  character(len=40)  :: bandnum, shortname
  character(len=200) :: longname
  character(len=40)  :: scheme_string
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  real :: PenSW_absorb_minthick ! A thickness that is used to absorb the remaining shortwave heat
                                ! flux when that flux drops below PEN_SW_FLUX_ABSORB [H ~> m or kg m-2]
  real :: PenSW_minthick_dflt ! The default for PenSW_absorb_minthick [m]
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags
  integer :: isd, ied, jsd, jed, nz, n
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, '')

! parameters for CHL_A routines
  call get_param(param_file, mdl, "VAR_PEN_SW", CS%var_pen_sw, &
                 "If true, use one of the CHL_A schemes specified by "//&
                 "OPACITY_SCHEME to determine the e-folding depth of "//&
                 "incoming short wave radiation.", default=.false.)

  CS%opacity_scheme = NO_SCHEME ; scheme_string = ''
  if (CS%var_pen_sw) then
    call get_param(param_file, mdl, "OPACITY_SCHEME", tmpstr, &
                 "This character string specifies how chlorophyll "//&
                 "concentrations are translated into opacities. Currently "//&
                 "valid options include:\n"//&
                 " \t\t  MANIZZA_05 - Use Manizza et al., GRL, 2005. \n"//&
                 " \t\t  MOREL_88 - Use Morel, JGR, 1988. \n"//&
                 " \t\t  OHLMANN_03 - Use Ohlmann, J Clim, 2003.", &
                 default=MANIZZA_05_STRING)
    if (len_trim(tmpstr) > 0) then
      tmpstr = uppercase(tmpstr)
      select case (tmpstr)
        case (MANIZZA_05_STRING)
          CS%opacity_scheme = MANIZZA_05 ; scheme_string = MANIZZA_05_STRING
        case (MOREL_88_STRING)
          CS%opacity_scheme = MOREL_88 ; scheme_string = MOREL_88_STRING
        case (OHLMANN_03_STRING)
          CS%opacity_scheme = OHLMANN_03 ; scheme_string = OHLMANN_03_STRING
        case default
          call MOM_error(FATAL, "opacity_init: #DEFINE OPACITY_SCHEME "//&
                                  trim(tmpstr) // "in input file is invalid.")
      end select
      call MOM_mesg('opacity_init: opacity scheme set to "'//trim(tmpstr)//'".', 5)
    endif
    if (CS%opacity_scheme == NO_SCHEME) then
      call MOM_error(WARNING, "opacity_init: No scheme has successfully "//&
               "been specified for the opacity.  Using the default MANIZZA_05.")
      CS%opacity_scheme = MANIZZA_05 ; scheme_string = MANIZZA_05_STRING
    endif

    call get_param(param_file, mdl, "BLUE_FRAC_SW", CS%blue_frac, &
                 "The fraction of the penetrating shortwave radiation "//&
                 "that is in the blue band.", default=0.5, units="nondim")
  else
    call get_param(param_file, mdl, "EXP_OPACITY_SCHEME", tmpstr, &
                 "This character string specifies which exponential "//&
                 "opacity scheme to utilize. Currently "//&
                 "valid options include:\n"//&
                 " \t\t  SINGLE_EXP - Single Exponent decay. \n"//&
                 " \t\t  DOUBLE_EXP - Double Exponent decay.", &
                 default=Single_Exp_String)!New default for "else" above (non-Chl scheme)
    if (len_trim(tmpstr) > 0) then
      tmpstr = uppercase(tmpstr)
      select case (tmpstr)
        case (SINGLE_EXP_STRING)
          CS%opacity_scheme = SINGLE_EXP ; scheme_string = SINGLE_EXP_STRING
        case (DOUBLE_EXP_STRING)
          CS%opacity_scheme = DOUBLE_EXP ; scheme_string = DOUBLE_EXP_STRING
      end select
      call MOM_mesg('opacity_init: opacity scheme set to "'//trim(tmpstr)//'".', 5)
    endif
    call get_param(param_file, mdl, "PEN_SW_SCALE", CS%pen_sw_scale, &
                 "The vertical absorption e-folding depth of the penetrating shortwave radiation.", &
                 units="m", default=0.0, scale=US%m_to_Z)
    !BGR/ Added for opacity_scheme==double_exp read in 2nd exp-decay and fraction
    if (CS%Opacity_scheme == DOUBLE_EXP ) then
      call get_param(param_file, mdl, "PEN_SW_SCALE_2ND", CS%pen_sw_scale_2nd, &
                 "The (2nd) vertical absorption e-folding depth of the "//&
                 "penetrating shortwave radiation (use if SW_EXP_MODE==double.)", &
                 units="m", default=0.0, scale=US%m_to_Z)
      call get_param(param_file, mdl, "SW_1ST_EXP_RATIO", CS%sw_1st_exp_ratio, &
                 "The fraction of 1st vertical absorption e-folding depth "//&
                 "penetrating shortwave radiation if SW_EXP_MODE==double.",&
                  units="nondim", default=0.0)
    elseif (CS%OPACITY_SCHEME == Single_Exp) then
      !/Else disable 2nd_exp scheme
      CS%pen_sw_scale_2nd = 0.0
      CS%sw_1st_exp_ratio = 1.0
    endif
    call get_param(param_file, mdl, "PEN_SW_FRAC", CS%pen_sw_frac, &
                 "The fraction of the shortwave radiation that penetrates "//&
                 "below the surface.", units="nondim", default=0.0)

  endif
  call get_param(param_file, mdl, "PEN_SW_NBANDS", optics%nbands, &
                 "The number of bands of penetrating shortwave radiation.", &
                 default=1)
  if (CS%Opacity_scheme == DOUBLE_EXP ) then
    if (optics%nbands /= 2) call MOM_error(FATAL, &
        "set_opacity: \Cannot use a double_exp opacity scheme with nbands!=2.")
  elseif (CS%Opacity_scheme == SINGLE_EXP ) then
    if (optics%nbands /= 1) call MOM_error(FATAL, &
        "set_opacity: \Cannot use a single_exp opacity scheme with nbands!=1.")
  elseif (CS%Opacity_scheme == OHLMANN_03 ) then
    if (optics%nbands /= 2) call MOM_error(FATAL, &
         "set_opacity: \OHLMANN_03 scheme requires nbands==2")
  endif

  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "OPTICS_ANSWER_DATE", optics%answer_date, &
                 "The vintage of the order of arithmetic and expressions in the optics calculations.  "//&
                 "Values below 20190101 recover the answers from the end of 2018, while "//&
                 "higher values use updated and more robust forms of the same expressions.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) optics%answer_date = max(optics%answer_date, 20230701)

  call get_param(param_file, mdl, "PEN_SW_FLUX_ABSORB", optics%PenSW_flux_absorb, &
                 "A minimum remaining shortwave heating rate that will be simply absorbed in "//&
                 "the next sufficiently thick layers for computational efficiency, instead of "//&
                 "continuing to penetrate.  The default, 2.5e-11 degC m s-1, is about 1e-4 W m-2 "//&
                 "or 0.08 degC m century-1, but 0 is also a valid value.", &
                 default=2.5e-11, units="degC m s-1", scale=US%degC_to_C*GV%m_to_H*US%T_to_s)

  if (optics%answer_date < 20190101) then ; PenSW_minthick_dflt = 0.001 ; else ; PenSW_minthick_dflt = 1.0 ; endif
  call get_param(param_file, mdl, "PEN_SW_ABSORB_MINTHICK", PenSW_absorb_minthick, &
                 "A thickness that is used to absorb the remaining penetrating shortwave heat "//&
                 "flux when it drops below PEN_SW_FLUX_ABSORB.", &
                 default=PenSW_minthick_dflt, units="m", scale=GV%m_to_H)
  optics%PenSW_absorb_Invlen = 1.0 / (PenSW_absorb_minthick + GV%H_subroundoff)

  if (.not.allocated(optics%min_wavelength_band)) &
    allocate(optics%min_wavelength_band(optics%nbands))
  if (.not.allocated(optics%max_wavelength_band)) &
    allocate(optics%max_wavelength_band(optics%nbands))

  if (CS%opacity_scheme == MANIZZA_05) then
    optics%min_wavelength_band(1) =0
    optics%max_wavelength_band(1) =550
    if (optics%nbands >= 2) then
      optics%min_wavelength_band(2)=550
      optics%max_wavelength_band(2)=700
    endif
    if (optics%nbands > 2) then
      do n=3,optics%nbands
        optics%min_wavelength_band(n) =700
        optics%max_wavelength_band(n) =2800
      enddo
    endif
  endif

  call get_param(param_file, mdl, "OPACITY_LAND_VALUE", CS%opacity_land_value, &
                 "The value to use for opacity over land. The default is "//&
                 "10 m-1 - a value for muddy water.", units="m-1", default=10.0, scale=US%Z_to_m)

  CS%warning_issued = .false.

  if (.not.allocated(optics%opacity_band)) &
    allocate(optics%opacity_band(optics%nbands,isd:ied,jsd:jed,nz), source=0.0)
  if (.not.allocated(optics%sw_pen_band)) &
    allocate(optics%sw_pen_band(optics%nbands,isd:ied,jsd:jed))
  allocate(CS%id_opacity(optics%nbands), source=-1)

  CS%id_sw_pen = register_diag_field('ocean_model', 'SW_pen', diag%axesT1, Time, &
      'Penetrating shortwave radiation flux into ocean', 'W m-2', conversion=US%QRZ_T_to_W_m2)
  CS%id_sw_vis_pen = register_diag_field('ocean_model', 'SW_vis_pen', diag%axesT1, Time, &
      'Visible penetrating shortwave radiation flux into ocean', 'W m-2', conversion=US%QRZ_T_to_W_m2)
  do n=1,optics%nbands
    write(bandnum,'(i3)') n
    shortname = 'opac_'//trim(adjustl(bandnum))
    longname = 'Opacity for shortwave radiation in band '//trim(adjustl(bandnum)) &
      // ', saved as L^-1 tanh(Opacity * L) for L = 10^-10 m'
    CS%id_opacity(n) = register_diag_field('ocean_model', shortname, diag%axesTL, Time, &
      longname, 'm-1', conversion=US%m_to_Z)
  enddo

  !! FOB
  if (CS%opacity_scheme == OHLMANN_03) then
     ! Set up the lookup table
     call init_ohlmann_table(optics)
  endif
  !! FOB

end subroutine opacity_init

!> Initialize the lookup table for Ohlmann solar penetration scheme.
!! Step size in Chl is a constant in log-space to make lookups easy.
!! Step size is fine enough that nearest neighbor lookup is sufficiently
!! accurate.
subroutine init_ohlmann_table(optics)

  implicit none

  type(optics_type), intent(inout) :: optics

  ! Local variables

  !! These are the data from Ohlmann (2003) Table 1a with additional
  !! values provided by C. Ohlmann and implemented in CESM-POP by B. Briegleb
  integer, parameter :: nval_tab1a = 31
  real, parameter, dimension(nval_tab1a) :: &
       chl_tab1a = (/                       &
       .001, .005, .01,  .02,               &
       .03,  .05,  .10,  .15,               &
       .20,  .25,  .30,  .35,               &
       .40,  .45,  .50,  .60,               &
       .70,  .80,  .90, 1.00,               &
       1.50, 2.00, 2.50, 3.00,              &
       4.00, 5.00, 6.00, 7.00,              &
       8.00, 9.00, 10.00  /)

  real, parameter, dimension(nval_tab1a) :: &
       a1_tab1a = (/                        &
       0.4421, 0.4451, 0.4488, 0.4563,      &
       0.4622, 0.4715, 0.4877, 0.4993,      &
       0.5084, 0.5159, 0.5223, 0.5278,      &
       0.5326, 0.5369, 0.5408, 0.5474,      &
       0.5529, 0.5576, 0.5615, 0.5649,      &
       0.5757, 0.5802, 0.5808, 0.5788,      &
       0.56965, 0.55638, 0.54091, 0.52442,  &
       0.50766, 0.49110, 0.47505  /)

  real, parameter, dimension(nval_tab1a) :: &
       a2_tab1a = (/                        &
       0.2981, 0.2963, 0.2940, 0.2894,      &
       0.2858, 0.2800, 0.2703, 0.2628,      &
       0.2571, 0.2523, 0.2481, 0.2444,      &
       0.2411, 0.2382, 0.2356, 0.2309,      &
       0.2269, 0.2235, 0.2206, 0.2181,      &
       0.2106, 0.2089, 0.2113, 0.2167,      &
       0.23357, 0.25504, 0.27829, 0.30274,  &
       0.32698, 0.35056, 0.37303 /)

  real, parameter, dimension(nval_tab1a) :: &
       b1_tab1a = (/                        &
       0.0287, 0.0301, 0.0319, 0.0355,      &
       0.0384, 0.0434, 0.0532, 0.0612,      &
       0.0681, 0.0743, 0.0800, 0.0853,      &
       0.0902, 0.0949, 0.0993, 0.1077,      &
       0.1154, 0.1227, 0.1294, 0.1359,      &
       0.1640, 0.1876, 0.2082, 0.2264,      &
       0.25808, 0.28498, 0.30844, 0.32932,  &
       0.34817, 0.36540, 0.38132 /)

  real, parameter, dimension(nval_tab1a) :: &
       b2_tab1a = (/                        &
       0.3192, 0.3243, 0.3306, 0.3433,      &
       0.3537, 0.3705, 0.4031, 0.4262,      &
       0.4456, 0.4621, 0.4763, 0.4889,      &
       0.4999, 0.5100, 0.5191, 0.5347,      &
       0.5477, 0.5588, 0.5682, 0.5764,      &
       0.6042, 0.6206, 0.6324, 0.6425,      &
       0.66172, 0.68144, 0.70086, 0.72144,  &
       0.74178, 0.76190, 0.78155 /)

  !! Make the table big enough so step size is smaller
  !! in log-space that any increment in Table 1a
  integer, parameter :: nval_lut=401
  real :: chl, log10chl_lut, w1, w2
  integer :: n,m,mm1,err

  allocate(optics%a1_lut(nval_lut),optics%b1_lut(nval_lut),&
       &   optics%a2_lut(nval_lut),optics%b2_lut(nval_lut),&
       &   stat=err)
  if ( err /= 0 ) then
     call MOM_error(FATAL,"init_ohlmann: Cannot allocate lookup table")
  endif

  optics%chl_min = chl_tab1a(1)
  optics%log10chl_min = log10(chl_tab1a(1))
  optics%log10chl_max = log10(chl_tab1a(nval_tab1a))
  optics%dlog10chl = (optics%log10chl_max - optics%log10chl_min)/(nval_lut-1)

  ! step through the lookup table
  m = 2
  do n=1,nval_lut
     log10chl_lut = optics%log10chl_min + (n-1)*optics%dlog10chl
     chl = 10.0**log10chl_lut
     chl = max(chl_tab1a(1),min(chl,chl_tab1a(nval_tab1a)))

     ! find interval in Table 1a (m-1,m]
     do while (chl > chl_tab1a(m))
        m = m + 1
     enddo
     mm1 = m-1

     ! interpolation weights
     w2 = (chl - chl_tab1a(mm1))/(chl_tab1a(m) - chl_tab1a(mm1))
     w1 = 1. - w2

     ! fill in the tables
     optics%a1_lut(n) = w1*a1_tab1a(mm1) + w2*a1_tab1a(m)
     optics%a2_lut(n) = w1*a2_tab1a(mm1) + w2*a2_tab1a(m)
     optics%b1_lut(n) = w1*b1_tab1a(mm1) + w2*b1_tab1a(m)
     optics%b2_lut(n) = w1*b2_tab1a(mm1) + w2*b2_tab1a(m)
  enddo

  return
end subroutine init_ohlmann_table

!> Get the partion of total solar into bands from Ohlmann lookup table
function lookup_ohlmann_swpen(chl,optics) result(A)

  implicit none

  real, intent(in) :: chl
  type(optics_type), intent(in) :: optics
  real, dimension(2) :: A

  ! Local variables

  real :: log10chl
  integer :: n

  ! Make sure we are in the table
  if (chl > optics%chl_min) then
    log10chl = min(log10(chl),optics%log10chl_max)
  else
    log10chl = optics%log10chl_min
  endif
  ! Do a nearest neighbor lookup
  n = nint( (log10chl - optics%log10chl_min)/optics%dlog10chl ) + 1

  A(1) = optics%a1_lut(n)
  A(2) = optics%a2_lut(n)

end function lookup_ohlmann_swpen

!> Get the opacity (decay scale) from Ohlmann lookup table
function lookup_ohlmann_opacity(chl,optics) result(B)

  implicit none
  real, intent(in) :: chl
  type(optics_type), intent(in) :: optics
  real, dimension(2) :: B

  ! Local variables
  real :: log10chl
  integer :: n

  ! Make sure we are in the table
  if (chl > optics%chl_min) then
    log10chl = min(log10(chl),optics%log10chl_max)
  else
    log10chl = optics%log10chl_min
  endif
  ! Do a nearest neighbor lookup
  n = nint( (log10chl - optics%log10chl_min)/optics%dlog10chl ) + 1

  B(1) = optics%b1_lut(n)
  B(2) = optics%b2_lut(n)

  return
end function lookup_ohlmann_opacity

subroutine opacity_end(CS, optics)
  type(opacity_CS)  :: CS     !< Opacity control structure
  type(optics_type) :: optics !< An optics type structure that should be deallocated.

  if (allocated(CS%id_opacity)) &
    deallocate(CS%id_opacity)
  if (allocated(optics%sw_pen_band)) &
    deallocate(optics%sw_pen_band)
  if (allocated(optics%opacity_band)) &
    deallocate(optics%opacity_band)
  if (allocated(optics%max_wavelength_band)) &
    deallocate(optics%max_wavelength_band)
  if (allocated(optics%min_wavelength_band)) &
       deallocate(optics%min_wavelength_band)
  if (allocated(optics%a1_lut)) deallocate(optics%a1_lut)
  if (allocated(optics%a2_lut)) deallocate(optics%a2_lut)
  if (allocated(optics%b1_lut)) deallocate(optics%b1_lut)
  if (allocated(optics%b2_lut)) deallocate(optics%b2_lut)
end subroutine opacity_end

!> \namespace mom_opacity
!!
!! opacity_from_chl:
!!   In this routine, the Morel (modified) or Manizza (modified)
!! schemes use the "blue" band in the parameterizations to determine
!! the e-folding depth of the incoming shortwave attenuation. The red
!! portion is lumped into the net heating at the surface.
!!
!! Morel, A., 1988: Optical modeling of the upper ocean in relation
!!   to its biogenous matter content (case-i waters)., J. Geo. Res.,
!!   93, 10,749-10,768.
!!
!! Manizza, M., C. LeQuere, A. J. Watson, and E. T. Buitenhuis, 2005:
!!  Bio-optical feedbacks among phytoplankton, upper ocean physics
!!  and sea-ice in a global model, Geophys. Res. Let., 32, L05603,
!!  doi:10.1029/2004GL020778.

!! Ohlmann, J.C., 2003: Ocean radiant heating in climate models.
!!   J. Climate, 16, 1337-1351, 2003.

end module MOM_opacity
