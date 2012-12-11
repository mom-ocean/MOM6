module adjustment_initialization

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : read_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, create_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use MOM_io, only : write_field, slasher, vardesc
use MOM_variables, only : thermo_var_ptrs, directories, ocean_OBC_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

implicit none ; private

#include <MOM_memory.h>

! -----------------------------------------------------------------------------
! Private (module-wise) parameters
! -----------------------------------------------------------------------------
integer, parameter :: IC_Z     = 0;     ! z coordinates
integer, parameter :: IC_RHO_L = 1;     ! layered isopycnals    
integer, parameter :: IC_RHO_C = 2;     ! continuous isopycnals
integer, parameter :: IC_SIGMA = 3;     ! sigma coordinates

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public adjustment_initialize_thickness
public adjustment_initialize_temperature_salinity 

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
! Initialization of thicknesses
!------------------------------------------------------------------------------
subroutine adjustment_initialize_thickness ( h, G, param_file )

  real, intent(out), dimension(NXMEM_,NYMEM_, NKMEM_) :: h
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file

! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!  This subroutine initializes the layer thicknesses to be uniform.
  real :: e0(SZK_(G)+1)   ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  real :: max_depth ! The maximum depths in m.
  integer :: i, j, k, is, ie, js, je, nz
  real    :: lenlon, lenlat
  real    :: x, y, yy, delta_S_strat, dSdz, delta_S, S_ref
  real    :: min_thickness, adjustment_width, adjustment_delta, adjustment_deltaS
  real    :: front_wave_amp, front_wave_length, front_wave_asym
  integer :: adjustment_ic
  real    :: target_values(SZK_(G)+1)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("MOM_initialization.F90, initialize_thickness_uniform: setting thickness")

  call read_param(param_file,"MAXIMUM_DEPTH",max_depth,.true.)
  call read_param (param_file,"LENLON",lenlon)
  call read_param (param_file,"LENLAT",lenlat)
  call read_param(param_file,"S_REF",S_ref,.true.)
  min_thickness = 1.0e-3; call read_param (param_file,"MIN_THICKNESS",min_thickness)
  adjustment_ic=-1; call read_param (param_file,"ADJUSTMENT_IC", adjustment_ic,.true.)
  adjustment_width=0.; call read_param (param_file,"ADJUSTMENT_WIDTH",adjustment_width,.true.)
  delta_S_strat = 0.0; call read_param(param_file,"DELTA_S_STRAT",delta_S_strat,.true.)
  adjustment_deltaS=0.; call read_param (param_file,"ADJUSTMENT_DELTAS",adjustment_deltaS,.true.)
  front_wave_amp=0.; call read_param (param_file,"FRONT_WAVE_AMP",front_wave_amp)
  front_wave_length=0.; call read_param (param_file,"FRONT_WAVE_LENGTH",front_wave_length)
  front_wave_asym=0.; call read_param (param_file,"FRONT_WAVE_ASYM",front_wave_asym)
 
  ! WARNING: this routine specifies the interface heights so that the last layer
  !          is vanished, even at maximum depth. In order to have a uniform
  !          layer distribution, use this line of code within the loop:
  !          e0(k) = -max_depth * real(k-1) / real(nz)
  !          To obtain a thickness distribution where the last layer is 
  !          vanished and the other thicknesses uniformly distributed, use:
  !          e0(k) = -max_depth * real(k-1) / real(nz-1)

  dSdz = -delta_S_strat/max_depth

  select case ( adjustment_ic )

    case ( IC_RHO_L, IC_RHO_C )
      if (delta_S_strat.ne.0.) then
        adjustment_delta = adjustment_deltaS / delta_S_strat * max_depth
        do k=1,nz+1
          e0(k) = adjustment_delta-(max_depth+2*adjustment_delta) * (real(k-1) / real(nz))
        enddo
      else
        adjustment_delta = 2.*max_depth
        do k=1,nz+1
          e0(k) = -(max_depth) * (real(k-1) / real(nz))
        enddo
      endif
      target_values(1)    = G%Rlay(1)+0.5*(G%Rlay(1)-G%Rlay(2))
      target_values(nz+1) = G%Rlay(nz)+0.5*(G%Rlay(nz)-G%Rlay(nz-1))
      do k = 2,nz
        target_values(k) = target_values(k-1) + ( G%Rlay(nz) - G%Rlay(1) ) / (nz-1)
      end do
      target_values = target_values - 1000.
      do j=js,je ; do i=is,ie
          if (front_wave_length.ne.0.) then
            y = ( 0.125 + G%geolath(i,j) / front_wave_length ) * ( 4. * acos(0.) )
            yy = 2. * ( G%geolath(i,j) - 0.5 * lenlat ) / adjustment_width
            yy = min(1.0, yy); yy = max(-1.0, yy)
            yy = yy * 2. * acos( 0. )
            y = front_wave_amp*sin(y) + front_wave_asym*sin(yy)
          else
            y = 0.
          endif
          x = ( ( G%geolonh(i,j) - 0.5 * lenlon ) + y ) / adjustment_width
          x = min(1.0, x); x = max(-1.0, x)
          x = x * acos( 0. )
          delta_S = adjustment_deltaS * 0.5 * (1. - sin( x ) )
          do k=2,nz
            if (dSdz.ne.0.) then
              eta1D(k) = ( target_values(k) - ( S_ref + delta_S ) ) / dSdz
            else
              eta1D(k) = e0(k) - (0.5*adjustment_delta) * sin( x )
            endif
            eta1D(k) = max( eta1D(k), -max_depth )
            eta1D(k) = min( eta1D(k), 0. )
          enddo
          eta1D(1)=0.; eta1D(nz+1)=-max_depth
          do k=nz,1,-1
            if (eta1D(k) > 0.) then
              eta1D(k) = max( eta1D(k+1) + min_thickness, 0. )
              h(i,j,k) = max( eta1D(k) - eta1D(k+1), min_thickness )
            elseif (eta1D(k) <= (eta1D(k+1) + min_thickness)) then
              eta1D(k) = eta1D(k+1) + min_thickness
              h(i,j,k) = min_thickness
            else
              h(i,j,k) = eta1D(k) - eta1D(k+1)
            endif
          enddo
      enddo ; enddo

    case ( IC_Z, IC_SIGMA )
      do k=1,nz+1
        eta1D(k) = -(max_depth) * (real(k-1) / real(nz))
        eta1D(k) = max(min(eta1D(k),0.),-max_depth)
      enddo
      do j=js,je ; do i=is,ie
        do k=nz,1,-1
            h(i,j,k) = eta1D(k) - eta1D(k+1)
        enddo
      enddo ; enddo

    case default
      call MOM_error(FATAL,"adjustment_initialize_thickness: "// &
      "Unrecognized i.c. setup - set ADJUSTMENT_IC")

  end select

end subroutine adjustment_initialize_thickness


!------------------------------------------------------------------------------
! Initialization of temperature and salinity
!------------------------------------------------------------------------------
subroutine adjustment_initialize_temperature_salinity ( T, S, h, G, param_file, &
                                                    eqn_of_state)
  real, dimension(NXMEM_,NYMEM_, NKMEM_), intent(out) :: T, S
  real, intent(in), dimension(NXMEM_,NYMEM_, NKMEM_)  :: h
  type(ocean_grid_type),               intent(in)  :: G
  type(param_file_type),               intent(in)  :: param_file
  type(EOS_type),                      pointer     :: eqn_of_state

  integer   :: i, j, k, is, ie, js, je, nz
  real      :: x, y, yy
  real      :: lenlon, lenlat
  real      :: max_depth
  integer   :: index_bay_z
  real      :: S_ref, T_ref         ! Reference salinity and temerature within
                                    ! surface layer
  real      :: S_range, T_range     ! Range of salinities and temperatures over the
                                    ! vertical
  real      :: xi0, xi1, dSdz, delta_S, delta_S_strat
  real      :: adjustment_width, adjustment_deltaS
  real       :: front_wave_amp, front_wave_length, front_wave_asym
  integer   :: adjustment_ic
  real      :: eta1d(SZK_(G)+1)
  
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call read_param(param_file,"LENLON",lenlon)
  call read_param(param_file,"LENLAT",lenlat)
  call read_param(param_file,"S_REF",S_ref,.true.)
  call read_param(param_file,"T_REF",T_ref,.true.)
  S_range = 2.0; call read_param(param_file,"S_RANGE",S_range,.false.)
  T_range = 0.0; call read_param(param_file,"T_RANGE",T_range,.false.)
  call read_param(param_file,"MAXIMUM_DEPTH",max_depth)
  call read_param(param_file,"ADJUSTMENT_IC",adjustment_ic)
  adjustment_width=0.; call read_param (param_file,"ADJUSTMENT_WIDTH",adjustment_width,.true.)
  adjustment_deltaS=0.; call read_param (param_file,"ADJUSTMENT_DELTAS",adjustment_deltaS,.true.)
  delta_S_strat = 0.0; call read_param(param_file,"DELTA_S_STRAT",delta_S_strat,.true.)
  front_wave_amp=0.; call read_param (param_file,"FRONT_WAVE_AMP",front_wave_amp)
  front_wave_length=0.; call read_param (param_file,"FRONT_WAVE_LENGTH",front_wave_length)
  front_wave_asym=0.; call read_param (param_file,"FRONT_WAVE_ASYM",front_wave_asym)

  dSdz = -delta_S_strat/max_depth
  T(:,:,:) = 0.0
  S(:,:,:) = 0.0
  
  ! Linear salinity profile
  select case ( adjustment_ic )

    case ( IC_Z, IC_SIGMA, IC_RHO_C )
      do j=js,je ; do i=is,ie
          eta1d(nz+1)=-G%D(i,j)
          do k=nz,1,-1
            eta1d(k)=eta1d(k+1)+h(i,j,k)
          enddo
          if (front_wave_length.ne.0.) then
            y = ( 0.125 + G%geolath(i,j) / front_wave_length ) * ( 4. * acos(0.) )
            yy = 2. * ( G%geolath(i,j) - 0.5 * lenlat ) / front_wave_length
            yy = min(1.0, yy); yy = max(-1.0, yy)
            yy = yy * 2. * acos( 0. )
            y = front_wave_amp*sin(y) + front_wave_asym*sin(yy)
          else
            y = 0.
          endif
          x = ( ( G%geolonh(i,j) - 0.5 * lenlon ) + y ) / adjustment_width
          x = min(1.0, x); x = max(-1.0, x)
          x = x * acos( 0. )
          delta_S = adjustment_deltaS * 0.5 * (1. - sin( x ) )
          do k=1,nz
            S(i,j,k) = S_ref + delta_S + 0.5 * ( eta1D(k)+eta1D(k+1) ) * dSdz
            x = abs(S(i,j,k) - 0.5*real(nz-1)/real(nz)*S_range)/S_range*real(2*nz)
            x = 1.-min(1., x)
            T(i,j,k) = x
         enddo
   !     x=sum(T(i,j,:)*h(i,j,:))
   !     T(i,j,:)=T(i,j,:)/x*(max_depth*1.5/real(nz))
      enddo ; enddo

    case ( IC_RHO_L ) 
      do k = 1,nz
        S(:,:,k) = S_ref + S_range * ( (real(k)-0.5) / real( nz ) )
   !    x = abs(S(1,1,k) - 0.5*real(nz-1)/real(nz)*S_range)/S_range*real(2*nz)
   !    x = 1.-min(1., x)
   !    T(:,:,k) = x
      end do 
    
    case default
      call MOM_error(FATAL,"adjustment_initialize_temperature_salinity: "// &
      "Unrecognized i.c. setup - set ADJUSTMENT_IC")

  end select

end subroutine adjustment_initialize_temperature_salinity

end module adjustment_initialization
