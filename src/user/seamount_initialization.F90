module seamount_initialization

use MOM_domains, only : sum_across_PEs
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : read_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, create_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use MOM_io, only : write_field, slasher, vardesc
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer, only : add_tracer_OBC_values, advect_tracer_CS
use MOM_variables, only : thermo_var_ptrs, directories, ocean_OBC_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

implicit none ; private

#include <MOM_memory.h>

! -----------------------------------------------------------------------------
! Private (module-wise) parameters
! -----------------------------------------------------------------------------
real, parameter :: seamount_delta = 0.5;
real, parameter :: seamount_length_scale = 20;  ! [km]

integer, parameter :: IC_Z     = 0; 
integer, parameter :: IC_RHO_L = 1;
integer, parameter :: IC_RHO_C = 2;
integer, parameter :: IC_SIGMA = 3;

integer, parameter :: RHO_LINEAR = 0;
integer, parameter :: RHO_EXP = 1;
integer, parameter :: RHO_PARABOLIC = 2;

integer, parameter :: seamount_ic = IC_SIGMA;
integer, parameter :: density_profile = RHO_LINEAR;


! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public seamount_initialize_topography
public seamount_initialize_thickness
public seamount_initialize_temperature_salinity 

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
! Initialization of topography
!------------------------------------------------------------------------------
subroutine seamount_initialize_topography ( D, G, param_file )
  ! Arguments 
  real, dimension(NIMEM_,NJMEM_), intent(out) :: D
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  
  ! Local variables 
  integer   :: i, j
  real      :: x;
  real      :: max_depth;
  real      :: delta;
  real      :: L;
  
  ! Maximum depth
  call read_param ( param_file, "MAXIMUM_DEPTH", max_depth );

  ! Domain extent in kilometers
  
  delta = seamount_delta;
  L = seamount_length_scale; 
  L = L / G%len_lon;
    
  do i=G%isc,G%iec 
    do j=G%jsc,G%jec 
    
      ! Compute normalized zonal coordinate (x=0 at center of domain)
      x = G%geoLonT(i,j) / G%len_lon - 0.5;
        
      D(i,j) = max_depth * ( 1.0 - delta * exp(-(x/L)**2) );

    enddo
  enddo

end subroutine seamount_initialize_topography


!------------------------------------------------------------------------------
! Initialization of thicknesses
!------------------------------------------------------------------------------
subroutine seamount_initialize_thickness ( h, G, param_file )

  real, intent(out), dimension(NIMEM_,NJMEM_, NKMEM_) :: h
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file

! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!  This subroutine initializes the layer thicknesses to be uniform.
  real :: e0(SZK_(G))     ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  real :: max_depth ! The maximum depths in m.
  integer :: i, j, k, is, ie, js, je, nz
  real    :: x;
  real    :: delta_h;
  real    :: min_thickness;

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("MOM_initialization.F90, initialize_thickness_uniform: setting thickness")

  call read_param(param_file,"MAXIMUM_DEPTH",max_depth,.true.)
  
  min_thickness = 1.0e-3; 
  call read_param ( param_file, "MIN_THICKNESS", min_thickness );
 
  ! WARNING: this routine specifies the interface heights so that the last layer
  !          is vanished, even at maximum depth. In order to have a uniform
  !          layer distribution, use this line of code within the loop:
  !          e0(k) = -max_depth * real(k-1) / real(nz)
  !          To obtain a thickness distribution where the last layer is 
  !          vanished and the other thicknesses uniformly distributed, use:
  !          e0(k) = -max_depth * real(k-1) / real(nz-1)
  do k=1,nz
    e0(k) = -max_depth * real(k-1) / real(nz)
  enddo

    
  select case ( seamount_ic )
    
  case ( IC_RHO_L )                 ! Initial thicknesses for isopycnal coordinates
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -1.0*G%bathyT(i,j)
      do k=nz,1,-1
        eta1D(k) = e0(k)
        if (eta1D(k) < (eta1D(k+1) + G%Angstrom_z)) then
          eta1D(k) = eta1D(k+1) + G%Angstrom_z
          h(i,j,k) = G%Angstrom_z
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
    enddo ; enddo

  case ( IC_Z, IC_RHO_C )                       ! Initial thicknesses for z coordinates
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -1.0*G%bathyT(i,j)
      do k=nz,1,-1
        eta1D(k) = e0(k)
        if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
          eta1D(k) = eta1D(k+1) + min_thickness
          h(i,j,k) = min_thickness
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
   enddo ; enddo

  case ( IC_SIGMA )             ! Initial thicknesses for sigma coordinates
    do j=js,je ; do i=is,ie
      delta_h = G%bathyT(i,j) / nz;
      h(i,j,:) = delta_h;
    end do ; end do 
end select
    

  call log_param(param_file, "seamount_uniform", "MAXIMUM_DEPTH",max_depth)

end subroutine seamount_initialize_thickness

!------------------------------------------------------------------------------
! Initialization of temperature and salinity
!------------------------------------------------------------------------------
subroutine seamount_initialize_temperature_salinity ( T, S, h, G, param_file, &
                                                      eqn_of_state)
                                                      
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: T, S
  real, intent(in), dimension(NIMEM_,NJMEM_, NKMEM_)  :: h
  type(ocean_grid_type),               intent(in)  :: G
  type(param_file_type),               intent(in)  :: param_file
  type(EOS_type),                      pointer     :: eqn_of_state

  integer   :: i, j, k, is, ie, js, je, nz
  real      :: max_depth
  real      :: xi0, xi1, dxi;
  real      :: r;
  
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  
  call read_param ( param_file, "MAXIMUM_DEPTH", max_depth );
 
  do j=js,je ; do i=is,ie
    
    xi0 = 0.0;
    do k = 1,nz
      xi1 = xi0 + h(i,j,k) / max_depth;
    
      if ( density_profile .eq. RHO_LINEAR ) then 
        ! ---------------------------
        ! Linear density profile
        ! ---------------------------
        S(i,j,k) = 34.0 + (xi0 + xi1);
      else if ( density_profile .eq. RHO_PARABOLIC ) then
        ! ---------------------------
        ! Parabolic density profile
        ! ---------------------------
        S(i,j,k) = 34.0 + (2.0 / 3.0) * (xi1**3 - xi0**3) / (xi1 - xi0);
      else if ( density_profile .eq. RHO_EXP ) then
        ! ---------------------------
        ! Exponential density profile
        ! ---------------------------
        r = 0.8;    ! small values give sharp profiles
        S(i,j,k) = 34.0 + r * (exp(xi1/r)-exp(xi0/r)) / (xi1 - xi0);
      end if    
      
      xi0 = xi1;

    enddo
  
  enddo ; enddo
  
end subroutine seamount_initialize_temperature_salinity

end module seamount_initialization
