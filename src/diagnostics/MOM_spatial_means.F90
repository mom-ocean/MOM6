!> Functions and routines to take area, volume, mass-weighted, layerwise, zonal or meridional means
module MOM_spatial_means

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms, only : EFP_type, operator(+), operator(-), assignment(=)
use MOM_coms, only : EFP_to_real, real_to_EFP, EFP_sum_across_PEs
use MOM_coms, only : reproducing_sum, reproducing_sum_EFP, EFP_to_real
use MOM_coms, only : query_EFP_overflow_error, reset_EFP_overflow_error
use MOM_coms, only : max_across_PEs, min_across_PEs
use MOM_error_handler, only : MOM_error, NOTE, WARNING, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public :: global_i_mean, global_j_mean
public :: global_area_mean, global_area_mean_u, global_area_mean_v, global_layer_mean
public :: global_area_integral
public :: global_volume_mean, global_mass_integral, global_mass_int_EFP
public :: adjust_area_mean_to_zero
public :: array_global_min_max

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.
! The functions in this module work with variables with arbitrary units, in which case the
! arbitrary rescaled units are indicated with [A ~> a], while the unscaled units are just [a].

contains

!> Return the global area mean of a variable, perhaps with a change of units. This uses reproducing sums.
function global_area_mean(var, G, scale, tmp_scale, unscale)
  type(ocean_grid_type),             intent(in)  :: G    !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)  :: var  !< The variable to average in arbitrary units [a],
                                                         !! or arbitrary rescaled units [A ~> a] if unscale
                                                         !! or tmp_scale is present
                                                         !! arbitrary, possibly rescaled units [A ~> a]
  real,                    optional, intent(in)  :: scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                         !! that converts it back to unscaled (e.g., mks)
                                                         !! units to enable the use of the reproducing sums
  real,                    optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the variable
                                                         !! that is reversed in the return value [a A-1 ~> 1],
                                                         !! or [b B-1 ~> 1] if unscale is also present.
  real,                    optional, intent(in)  :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                         !! that converts it back to unscaled (e.g., mks)
                                                         !! units to enable the use of the reproducing sums, or
                                                         !! a factor converting between rescaled units if
                                                         !! tmp_scale is also present [B A-1 ~> b a-1].
                                                         !! Here scale and unscale are synonymous, but unscale
                                                         !! is preferred and takes precedence.
  real :: global_area_mean  ! The mean of the variable in arbitrary unscaled units [a] or scaled units [A ~> a]
                            ! or [B ~> b], depending on which optional arguments are provided

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums.
  ! [A ~> a] and [B ~> b] are the same units unless tmp_scale and unscale are both present.
  real :: tmpForSumming(SZI_(G),SZJ_(G))  ! An unscaled cell integral in [a L2 ~> a m2] or a
                                          ! scaled cell integral in [A L2 ~> a m2] or [B L2 ~> b m2]
  real :: scalefac   ! A scaling factor for the variable that is not reversed [a A-1 ~> 1] or [B A-1 ~> b a-1] or [1]
  real :: temp_scale ! A temporary scaling factor [a A-1 ~> 1] or [b B-1 ~> 1] or [1]
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  temp_scale = 1.0 ; if (present(tmp_scale)) temp_scale = tmp_scale

  scalefac = 1.0
  if (present(unscale)) then ; scalefac = unscale
  elseif (present(scale)) then ; scalefac = scale ; endif

  tmpForSumming(:,:) = 0.
  do j=js,je ; do i=is,ie
    tmpForSumming(i,j) = var(i,j) * (scalefac * G%areaT(i,j) * G%mask2dT(i,j))
  enddo ; enddo

  global_area_mean = reproducing_sum(tmpForSumming, unscale=temp_scale*G%US%L_to_m**2) * G%IareaT_global

end function global_area_mean

!> Return the global area mean of a variable. This uses reproducing sums.
function global_area_mean_v(var, G, tmp_scale)
  type(ocean_grid_type),             intent(in)  :: G    !< The ocean's grid structure
  real, dimension(SZI_(G),SZJB_(G)), intent(in)  :: var  !< The variable to average in arbitrary
                                                         !! units [a], or arbitrary rescaled units
                                                         !! [A ~> a] if tmp_scale is present
  real,                    optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the
                                                         !! variable that converts it back to unscaled
                                                         !! (e.g., mks) units to enable the use of the
                                                         !! reproducing sums [a A-1 ~> 1], but is reversed
                                                         !! before output so that the return value has
                                                         !! the same units as var

  real :: global_area_mean_v  ! The mean of the variable in the same arbitrary units as var [A ~> a]

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  ! [A ~> a] and [B ~> b] are the same unless tmp_scale and unscale are both present.
  real, dimension(SZI_(G),SZJ_(G)) :: tmpForSumming ! An unscaled cell integral [A L2 ~> a m2]
  real :: temp_scale ! A temporary scaling factor [a A-1 ~> 1] or [1]
  integer :: i, j, is, ie, js, je, isB, ieB, jsB, jeB

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isB = G%iscB ; ieB = G%iecB ; jsB = G%jscB ; jeB = G%jecB

  temp_scale = 1.0 ; if (present(tmp_scale)) temp_scale = tmp_scale

  tmpForSumming(:,:) = 0.
  do j=js,je ; do i=is,ie
    tmpForSumming(i,j) = G%areaT(i,j) * &
             (var(i,J) * G%mask2dCv(i,J) + var(i,J-1) * G%mask2dCv(i,J-1)) / &
             max(1.e-20, G%mask2dCv(i,J)+G%mask2dCv(i,J-1))
  enddo ; enddo
  global_area_mean_v = reproducing_sum(tmpForSumming, unscale=G%US%L_to_m**2*temp_scale) * G%IareaT_global

end function global_area_mean_v

!> Return the global area mean of a variable on U grid. This uses reproducing sums.
function global_area_mean_u(var, G, tmp_scale)
  type(ocean_grid_type),             intent(in)  :: G    !< The ocean's grid structure
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)  :: var  !< The variable to average in arbitrary
                                                         !! units [a], or arbitrary rescaled units
                                                         !! [A ~> a] if tmp_scale is present
  real,                    optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the
                                                         !! variable that converts it back to unscaled
                                                         !! (e.g., mks) units to enable the use of the
                                                         !! reproducing sums [a A-1 ~> 1], but is reversed
                                                         !! before output so that the return value has
                                                         !! the same units as var
  real :: global_area_mean_u  ! The mean of the variable in the same arbitrary units as var [A ~> a]

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  real, dimension(SZI_(G),SZJ_(G)) :: tmpForSumming ! An unscaled cell integral [A L2 ~> a m2]
  real :: temp_scale ! A temporary scaling factor [a A-1 ~> 1] or [1]
  integer :: i, j, is, ie, js, je, isB, ieB, jsB, jeB

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isB = G%iscB ; ieB = G%iecB ; jsB = G%jscB ; jeB = G%jecB

  temp_scale = 1.0 ; if (present(tmp_scale)) temp_scale = tmp_scale

  tmpForSumming(:,:) = 0.
  do j=js,je ; do i=is,ie
    tmpForSumming(i,j) = G%areaT(i,j) * &
             (var(I,j) * G%mask2dCu(I,j) + var(I-1,j) * G%mask2dCu(I-1,j)) / &
             max(1.e-20, G%mask2dCu(I,j)+G%mask2dCu(I-1,j))
  enddo ; enddo
  global_area_mean_u = reproducing_sum(tmpForSumming, unscale=G%US%L_to_m**2*temp_scale) * G%IareaT_global

end function global_area_mean_u

!> Return the global area integral of a variable using reproducing sums, perhaps with a change of units.
!! By default the integral uses the masked area from the grid, but an alternate could be used instead.
!! The presence of the optional tmp_scale argument determines whether the returned value is in scaled
!! (if it is present) or unscaled units for both the variable itself and for the area in the integral.
function global_area_integral(var, G, scale, area, tmp_scale, unscale)
  type(ocean_grid_type),            intent(in)  :: G     !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: var   !< The variable to integrate in arbitrary units [a],
                                                         !! or arbitrary rescaled units [A ~> a] if unscale
                                                         !! or tmp_scale is present
  real,                   optional, intent(in)  :: scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                         !! that converts it back to unscaled (e.g., mks)
                                                         !! units to enable the use of the reproducing sums
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: area !< The alternate area to use, including
                                                         !! any required masking [L2 ~> m2].
  real,                   optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the variable
                                                         !! that is reversed in the return value [a A-1 ~> 1],
                                                         !! or [b B-1 ~> 1] if unscale is also present.
  real,                   optional, intent(in)  :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                         !! that converts it back to unscaled (e.g., mks)
                                                         !! units to enable the use of the reproducing sums, or
                                                         !! a factor converting between rescaled units if
                                                         !! tmp_scale is also present [B A-1 ~> b a-1].
                                                         !! Here scale and unscale are synonymous, but unscale
                                                         !! is preferred and takes precedence if both are present.
  real :: global_area_integral !< The returned area integral, usually in the units of var times an area,
                               !! [a m2] or [A L2 ~> a m2] or [B L2 ~> b m2], depending on which optional
                               !! arguments are provided

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums.
  ! [A ~> a] and [B ~> b] are the same units unless tmp_scale and unscale are both present.
  real, dimension(SZI_(G),SZJ_(G)) :: tmpForSumming ! An unscaled cell integral in [a m2] or
                    ! a scaled cell integral in [B L2 ~> b m2] or other units as indicated below
  real :: scalefac  ! An overall scaling factor for the areas and variable, in units of [a m2 A-1 L-2 ~> 1]
                    ! or [1] or [B m2 A-1 L-2 ~> b a-1] or [B A-1 ~> b a-1] depending on which
                    ! optional arguments are present.
  !_______________________________________________________________________________________________
  ! Table of units of scalefac and tmpForSumming, depending on the presence of optional arguments |
  !_______________________________________________________________________________________________|
  ! present(tmp_scale) | present(unscale) | scalefac units          | tmpForSumming units         |
  !____________________|__________________|_________________________|_____________________________!
  !      True          |      True        | [B A-1 ~> b a-1]        | [B L2 ~> b m2]              |
  !      True          |      False       | [1]                     | [A L2 ~> a m2]              |
  !      False         |      True        | [a m2 A-1 L-2 ~> b a-1] | [a m2]                      |
  !      False         |      False       | [m2 L-2 ~> 1]           | [a m2]                      |
  !____________________|__________________|_________________________|_____________________________!
  real :: temp_scale ! A temporary scaling factor [a m2 L-2 A-1 ~> 1] or [b m2 L-2 B-1 ~> 1] or [1]
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (present(tmp_scale)) then
    temp_scale = G%US%L_to_m**2 * tmp_scale  ! Units of [a m2 A-1 L-2 ~> 1] or [b m2 B-1 L-2 ~> 1]
    scalefac = 1.0
  else
    temp_scale = 1.0
    scalefac = G%US%L_to_m**2
  endif
  if (present(unscale)) then ; scalefac = scalefac * unscale
  elseif (present(scale)) then ; scalefac = scalefac * scale ; endif

  tmpForSumming(:,:) = 0.
  if (present(area)) then
    do j=js,je ; do i=is,ie
      tmpForSumming(i,j) = var(i,j) * (scalefac * area(i,j))
    enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      tmpForSumming(i,j) = var(i,j) * (scalefac * G%areaT(i,j) * G%mask2dT(i,j))
    enddo ; enddo
  endif

  global_area_integral = reproducing_sum(tmpForSumming, unscale=temp_scale)

end function global_area_integral

!> Return the layerwise global thickness-weighted mean of a variable. This uses reproducing sums.
function global_layer_mean(var, h, G, GV, scale, tmp_scale, unscale)
  type(ocean_grid_type),                     intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: var  !< The variable to average in arbitrary units [a],
                                                                 !! or arbitrary rescaled units [A ~> a] if unscale
                                                                 !! or tmp_scale is present
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real,                            optional, intent(in)  :: scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                                 !! that converts it back to unscaled (e.g., mks)
                                                                 !! units to enable the use of the reproducing sums
  real,                            optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the
                                                                 !! variable for use in the reproducing sums
                                                                 !! that is reversed in the return value [a A-1 ~> 1],
                                                                 !! or [b B-1 ~> 1] if unscale is also present.
  real,                            optional, intent(in)  :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                                 !! that converts it back to unscaled (e.g., mks)
                                                                 !! units to enable the use of the reproducing sums, or
                                                                 !! a factor converting between rescaled units if
                                                                 !! tmp_scale is also present [B A-1 ~> b a-1].
                                                                 !! Here scale and unscale are synonymous, but unscale
                                                                 !! is preferred and takes precedence.
  real, dimension(SZK_(GV)) :: global_layer_mean  !< The mean of the variable in the arbitrary scaled [A ~> a]
                                                  !! or [B ~> b] or unscaled [a] units of var, depending on which
                                                  !! optional arguments are provided

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  real :: tmpForSumming(G%isc:G%iec,G%jsc:G%jec,SZK_(GV)) ! An unscaled cell integral in [L2 a m ~> a m3] or
                                                    ! [L2 a kg m-2 ~> a kg] or a scaled cell integral in
                                                    ! [L2 B m ~> b m3] or [L2 B m ~> b m3] or other units
                                                    ! as indicated the table below.
  real :: weight(G%isc:G%iec,G%jsc:G%jec,SZK_(GV))  ! The volume or mass of each cell, depending on whether
                                                    ! the model is Boussinesq, used as a weight [L2 m ~> m3]
                                                    ! or [L2 kg m-2 ~> kg]
  real :: scalefac   ! A scaling factor for the variable [a A-1 ~> 1] or [B A-1 ~> b a-1] or [1]
  !__________________________________________________________________________________________________
  ! Units of weight, scalefac and tmpForSumming, depending on the presence of optional arguments    |
  !_________________________________________________________________________________________________|
  ! Boussinesq | tmp_scale | unscale |  weight units     | scalefac units   | tmpForSumming units   |
  !            | present   | present |                   |                  |                       |
  !____________|___________|_________|___________________|__________________|_______________________!
  !   True     |  True     |  True   | [L2 m ~> m3]      | [B A-1 ~> b a-1] | [B L2 m ~> b m3]      |
  !   True     |  True     |  False  | [L2 m ~> m3]      | [1]              | [A L2 m ~> a m3]      |
  !   True     |  False    |  True   | [L2 m ~> m3]      | [a A-1 ~> 1]     | [L2 a m ~> a m3]      |
  !   True     |  False    |  False  | [L2 m ~> m3]      | [1]              | [L2 a m ~> a m3]      |
  !   False    |  True     |  True   | [L2 kg m-2 ~> kg] | [B A-1 ~> b a-1] | [B L2 kg m-2 ~> b kg] |
  !   False    |  True     |  False  | [L2 kg m-2 ~> kg] | [1]              | [A L2 kg m-2 ~> a kg] |
  !   False    |  False    |  True   | [L2 kg m-2 ~> kg] | [a A-1 ~> 1]     | [L2 a kg m-2 ~> a kg] |
  !   False    |  False    |  False  | [L2 kg m-2 ~> kg] | [1]              | [L2 a kg m-2 ~> a kg] |
  !____________|___________|_________|___________________|__________________|_______________________!
  type(EFP_type) :: laysums(2*SZK_(GV)) ! A vector of sums with heterogeneous meanings, with the first
                                        ! half being the tracer integrals in [b m3] or [b kg] and the
                                        ! second half being the summed weights in [m3] or [kg]
  real :: global_temp_scalar   ! The global integral of the tracer over all
                               ! layers [L2 a m ~> a m3] or [L2 a kg m-2 ~> a kg]
  real :: global_weight_scalar ! The global integral of the volume or mass over all
                               ! layers [L2 m ~> m3] or [L2 kg m-2 ~> kg]
  real :: temp_scale ! A temporary scaling factor [a A-1 ~> 1] or [b B-1 ~> 1] or [1]
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  temp_scale = 1.0 ; if (present(tmp_scale)) temp_scale = tmp_scale

  scalefac = 1.0
  if (present(unscale)) then ; scalefac = unscale
  elseif (present(scale)) then ; scalefac = scale ; endif
  tmpForSumming(:,:,:) = 0. ; weight(:,:,:) = 0.

  do k=1,nz ; do j=js,je ; do i=is,ie
    weight(i,j,k)  =  (GV%H_to_MKS * h(i,j,k)) * (G%areaT(i,j) * G%mask2dT(i,j))
    tmpForSumming(i,j,k) =  scalefac * var(i,j,k) * weight(i,j,k)
  enddo ; enddo ; enddo

  global_temp_scalar = reproducing_sum(tmpForSumming, EFP_lay_sums=laysums(1:nz), only_on_PE=.true., &
                                       unscale=temp_scale*G%US%L_to_m**2)
  global_weight_scalar = reproducing_sum(weight, EFP_lay_sums=laysums(nz+1:2*nz), only_on_PE=.true., &
                                         unscale=G%US%L_to_m**2)
  call EFP_sum_across_PEs(laysums, 2*nz)

  ! Note that temp_scale appears in the denominator here because the variables returned via the
  ! EFP_lay_sums arguments to reproducing sums stay in unscaled mks units.
  do k=1,nz
    global_layer_mean(k) = EFP_to_real(laysums(k)) / (temp_scale*EFP_to_real(laysums(nz+k)))
  enddo

end function global_layer_mean

!> Find the global thickness-weighted mean of a variable. This uses reproducing sums.
function global_volume_mean(var, h, G, GV, scale, tmp_scale, unscale)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: var  !< The variable to average in arbitrary units [a],
                                               !! or arbitrary rescaled units [A ~> a] if unscale
                                               !! or tmp_scale is present
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real,          optional, intent(in)  :: scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                               !! that converts it back to unscaled (e.g., mks)
                                               !! units to enable the use of the reproducing sums
  real,          optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the
                                               !! variable that is reversed in the return value [a A-1 ~> 1],
                                               !! or [b B-1 ~> 1] if unscale is also present.
  real,          optional, intent(in)  :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                               !! that converts it back to unscaled (e.g., mks)
                                               !! units to enable the use of the reproducing sums, or
                                               !! a factor converting between rescaled units if
                                               !! tmp_scale is also present [B A-1 ~> b a-1].
                                               !! Here scale and unscale are synonymous, but unscale
                                               !! is preferred and takes precedence if both are present.
  real :: global_volume_mean  !< The thickness-weighted average of var in the arbitrary scaled [A ~> a] or [B ~> b] or
                              !! unscaled [a] units of var, depending on which optional arguments are provided

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  ! [A ~> a] and [B ~> b] are the same units unless tmp_scale and unscale are both present.
  real :: temp_scale ! A temporary scaling factor [a A-1 ~> 1] or [b B-1 ~> 1] or [1]
  real :: scalefac   ! A scaling factor for the variable [a A-1 ~> 1] or [B A-1 ~> b a-1] or [1]
  real :: weight_here ! The volume or mass of a grid cell [L2 m ~> m3] or [L2 kg m-2 ~> kg]
  real, dimension(SZI_(G),SZJ_(G)) :: tmpForSumming ! The volume or mass integral of the variable in a column
                                                  ! [B L2 m ~> b m3] or [B L2 kg m-2 ~> b kg] or
                                                  ! [L2 a m ~> a m3] or [L2 a kg m-2 ~> a kg]
  real, dimension(SZI_(G),SZJ_(G)) :: sum_weight  ! The volume or mass of each column of water
                                                  ! [L2 m ~> m3] or [L2 kg m-2 ~> kg]
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  temp_scale = 1.0 ; if (present(tmp_scale)) temp_scale = tmp_scale

  scalefac = 1.0
  if (present(unscale)) then ; scalefac = unscale
  elseif (present(scale)) then ; scalefac = scale ; endif
  tmpForSumming(:,:) = 0. ; sum_weight(:,:) = 0.

  do k=1,nz ; do j=js,je ; do i=is,ie
    weight_here  =  (GV%H_to_MKS * h(i,j,k)) * (G%areaT(i,j) * G%mask2dT(i,j))
    tmpForSumming(i,j) = tmpForSumming(i,j) + scalefac * var(i,j,k) * weight_here
    sum_weight(i,j) = sum_weight(i,j) + weight_here
  enddo ; enddo ; enddo
  global_volume_mean = (reproducing_sum(tmpForSumming, unscale=temp_scale*G%US%L_to_m**2)) / &
                       (reproducing_sum(sum_weight, unscale=G%US%L_to_m**2))

end function global_volume_mean


!> Find the global mass-weighted integral of a variable.  The presence of the optional tmp_scale
!! argument determines whether the returned value is in scaled (if it is present) or unscaled units
!! for both the variable itself and for the mass in the integral.  This function uses reproducing sums.
function global_mass_integral(h, G, GV, var, on_PE_only, scale, tmp_scale, unscale)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)  :: var  !< The variable to integrate in arbitrary units [a],
                                               !! or arbitrary rescaled units [A ~> a] if unscale
                                               !! or tmp_scale is present
  logical,       optional, intent(in)  :: on_PE_only  !< If present and true, the sum is only done
                                               !! on the local PE, and it is _not_ order invariant.
  real,          optional, intent(in)  :: scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                               !! that converts it back to unscaled (e.g., mks)
                                               !! units to enable the use of the reproducing sums
  real,          optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the variable
                                               !! that is reversed in the return value [a A-1 ~> 1],
                                               !! or [b B-1 ~> 1] if unscale is also present.
  real,          optional, intent(in)  :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                               !! that converts it back to unscaled (e.g., mks)
                                               !! units to enable the use of the reproducing sums, or
                                               !! a factor converting between rescaled units if
                                               !! tmp_scale is also present [B A-1 ~> b a-1].
                                               !! Here scale and unscale are synonymous, but unscale
                                               !! is preferred and takes precedence if both are present.
  real :: global_mass_integral  !< The mass-weighted integral of var (or 1) in kg times the arbitrary
                                !! units of var [kg a] or in [R Z L2 A ~> kg a] if tmp_scale is present
                                !! or [R Z L2 B ~> kg b] if both unscale and tmp_scale are present

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  ! [A ~> a] and [B ~> b] are the same units unless tmp_scale and unscale are both present.
  real :: tmpForSumming(SZI_(G),SZJ_(G)) ! The mass-weighted integral of the variable in a column in
                     ! [kg a] or [kg] or if tmp_scale is present in [B R Z L2 ~> kg b] or
                     ! [A R Z L2 !> kg m] or [R Z L2 ~> kg]
  real :: scalefac   ! An overall scaling factor for the cell mass and variable in [a kg A-1 R-1 Z-1 L-2 ~> 1]
                     ! or [kg R-1 Z-1 L-2 ~> 1] or [1] or [B A-1 ~> b a-1] if tmp_scale is present.
  real :: temp_scale ! A temporary scaling factor [1] or if tmp_scale is present this could be in
                     ! [kg a R-1 Z-1 L-2 A-1 ~> 1] or [kg b R-1 Z-1 L-2 B-1 ~> 1] or [kg R-1 Z-1 L-2 ~> 1]
  !_______________________________________________________________________________________
  ! Units of scalefac and tmpForSumming, depending on the presence of optional arguments |
  !______________________________________________________________________________________|
  ! var     | tmp_scale | unscale |   scalefac units            | tmpForSumming units    |
  ! present | present   | present |                             |                        |
  !_________|___________|_________|_____________________________|________________________!
  !  True   |  True     |  True   | [B A-1 ~> b a-1]            | [B R Z L2 ~> b kg]     |
  !  True   |  True     |  False  | [1]                         | [A R Z L2 ~> a kg]     |
  !  True   |  False    |  True   | [a kg A-1 R-1 Z-1 L-2 ~> 1] | [a kg]                 |
  !  True   |  False    |  False  | [kg R-1 Z-1 L-2 ~> 1]       | [a kg]                 |
  !  False  |  True     |  either | [1]                         | [R Z L2 ~> kg]         |
  !  False  |  False    |  either | [kg R-1 Z-1 L-2 ~> 1]       | [kg]                   |
  !_________|___________|_________|_____________________________|________________________!
  logical :: global_sum ! If true do the sum globally, but if false only do the sum on the current PE.
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (present(tmp_scale)) then
    temp_scale = G%US%RZL2_to_kg * tmp_scale
    if (.not.present(var)) temp_scale = G%US%RZL2_to_kg
    scalefac = 1.0
  else
    temp_scale = 1.0
    scalefac = G%US%RZL2_to_kg
  endif
  if (present(var)) then
    if (present(unscale)) then ; scalefac = scalefac * unscale
    elseif (present(scale)) then ; scalefac = scalefac * scale ; endif
  endif

  tmpForSumming(:,:) = 0.0
  if (present(var)) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      tmpForSumming(i,j) = tmpForSumming(i,j) + var(i,j,k) * &
                ((GV%H_to_RZ * h(i,j,k)) * (scalefac*G%areaT(i,j) * G%mask2dT(i,j)))
    enddo ; enddo ; enddo
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      tmpForSumming(i,j) = tmpForSumming(i,j) + &
                ((GV%H_to_RZ * h(i,j,k)) * (scalefac*G%areaT(i,j) * G%mask2dT(i,j)))
    enddo ; enddo ; enddo
  endif
  global_sum = .true. ; if (present(on_PE_only)) global_sum = .not.on_PE_only
  if (global_sum) then
    global_mass_integral = reproducing_sum(tmpForSumming, unscale=temp_scale)
  else
    global_mass_integral = 0.0
    do j=js,je ; do i=is,ie
      global_mass_integral = global_mass_integral + tmpForSumming(i,j)
    enddo ; enddo
  endif

end function global_mass_integral

!> Find the global mass-weighted order invariant integral of a variable in mks units,
!! returning the value as an EFP_type. This uses reproducing sums.
function global_mass_int_EFP(h, G, GV, var, on_PE_only, scale, unscale)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)  :: var  !< The variable to integrate in arbitrary units [a],
                                               !! or arbitrary rescaled units [A ~> a] if unscale
                                               !! is present
  logical,       optional, intent(in)  :: on_PE_only  !< If present and true, the sum is only done
                                               !! on the local PE, but it is still order invariant.
  real,          optional, intent(in)  :: scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                               !! that converts it back to unscaled (e.g., mks)
                                               !! units to enable the use of the reproducing sums
  real,          optional, intent(in)  :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                               !! that converts it back to unscaled (e.g., mks)
                                               !! units to enable the use of the reproducing sums.
                                               !! Here scale and unscale are synonymous, but unscale
                                               !! is preferred and takes precedence if both are present.
  type(EFP_type) :: global_mass_int_EFP  !< The mass-weighted integral of var (or 1) in
                                         !! kg times the arbitrary units of var [kg a]

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  real :: tmpForSum(SZI_(G),SZJ_(G)) ! The mass-weighted integral of the variable in a column [kg a] or [kg]
  real :: scalefac  ! An overall scaling factor for the cell mass and variable [a kg A-1 H-1 L-2 ~> kg m-3 or 1]
  integer :: i, j, k, is, ie, js, je, nz, isr, ier, jsr, jer

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isr = is - (G%isd-1) ; ier = ie - (G%isd-1) ; jsr = js - (G%jsd-1) ; jer = je - (G%jsd-1)

  scalefac = GV%H_to_kg_m2 * G%US%L_to_m**2
  if (present(unscale)) then ; scalefac = unscale * scalefac
  elseif (present(scale)) then ; scalefac = scale * scalefac ; endif

  tmpForSum(:,:) = 0.0
  if (present(var)) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      tmpForSum(i,j) = tmpForSum(i,j) + var(i,j,k) * &
                ((scalefac * h(i,j,k)) * (G%areaT(i,j) * G%mask2dT(i,j)))
    enddo ; enddo ; enddo
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      tmpForSum(i,j) = tmpForSum(i,j) + &
                ((scalefac * h(i,j,k)) * (G%areaT(i,j) * G%mask2dT(i,j)))
    enddo ; enddo ; enddo
  endif

  global_mass_int_EFP = reproducing_sum_EFP(tmpForSum, isr, ier, jsr, jer, only_on_PE=on_PE_only)

end function global_mass_int_EFP


!> Determine the global mean of a field along rows of constant i, returning it
!! in a 1-d array using the local indexing. This uses reproducing sums.
subroutine global_i_mean(array, i_mean, G, mask, scale, tmp_scale, unscale)
  type(ocean_grid_type),            intent(inout) :: G     !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: array !< The variable to integrate in arbitrary units [a],
                                                           !! or arbitrary rescaled units [A ~> a] if unscale
                                                           !! is present
  real, dimension(SZJ_(G)),         intent(out)   :: i_mean !< Global mean of array along its i-axis [a] or [A ~> a]
  real, dimension(SZI_(G),SZJ_(G)), &
                          optional, intent(in)    :: mask  !< An array used for weighting the i-mean [nondim]
  real,                   optional, intent(in)    :: scale !< A rescaling factor for the output variable [a A-1 ~> 1]
                                                           !! that converts it back to unscaled (e.g., mks)
                                                           !! units to enable the use of the reproducing sums
  real,                   optional, intent(in)    :: tmp_scale !< A rescaling factor for the internal
                                                           !! calculations that is removed from the output [a A-1 ~> 1]
  real,                   optional, intent(in)    :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                           !! that converts it back to unscaled (e.g., mks)
                                                           !! units to enable the use of the reproducing sums.
                                                           !! Here scale and unscale are synonymous, but unscale
                                                           !! is preferred and takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  type(EFP_type), allocatable, dimension(:) :: asum      ! The masked sum of the variable in each row [a]
  type(EFP_type), allocatable, dimension(:) :: mask_sum  ! The sum of the mask values in each row [nondim]
  real :: scalefac   ! A scaling factor for the variable [a A-1 ~> 1]
  real :: rescale    ! A factor for redoing any internal rescaling before output [A a-1 ~> 1]
  real :: mask_sum_r ! The sum of the mask values in a row [nondim]
  integer :: is, ie, js, je, idg_off, jdg_off
  integer :: i, j

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  idg_off = G%idg_offset ; jdg_off = G%jdg_offset

  scalefac = 1.0
  if (present(unscale)) then ; scalefac = unscale
  elseif (present(scale)) then ; scalefac = scale ; endif

  rescale = 1.0
  if (present(tmp_scale)) then ; if (tmp_scale /= 0.0) then
    scalefac = scalefac * tmp_scale ; rescale = 1.0 / tmp_scale
  endif ; endif
  call reset_EFP_overflow_error()

  allocate(asum(G%jsg:G%jeg))
  if (present(mask)) then
    allocate(mask_sum(G%jsg:G%jeg))

    do j=G%jsg,G%jeg
      asum(j) = real_to_EFP(0.0) ; mask_sum(j) = real_to_EFP(0.0)
    enddo

    do j=js,je ; do i=is,ie
      asum(j+jdg_off) = asum(j+jdg_off) + real_to_EFP(scalefac*array(i,j)*mask(i,j))
      mask_sum(j+jdg_off) = mask_sum(j+jdg_off) + real_to_EFP(mask(i,j))
    enddo ; enddo

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_i_mean overflow error occurred before sums across PEs.")

    call EFP_sum_across_PEs(asum(G%jsg:G%jeg), G%jeg-G%jsg+1)
    call EFP_sum_across_PEs(mask_sum(G%jsg:G%jeg), G%jeg-G%jsg+1)

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_i_mean overflow error occurred during sums across PEs.")

    do j=js,je
      mask_sum_r = EFP_to_real(mask_sum(j+jdg_off))
      if (mask_sum_r == 0.0 ) then ; i_mean(j) = 0.0 ; else
        i_mean(j) = EFP_to_real(asum(j+jdg_off)) / mask_sum_r
      endif
    enddo

    deallocate(mask_sum)
  else
    do j=G%jsg,G%jeg ; asum(j) = real_to_EFP(0.0) ; enddo

    do j=js,je ; do i=is,ie
      asum(j+jdg_off) = asum(j+jdg_off) + real_to_EFP(scalefac*array(i,j))
    enddo ; enddo

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_i_mean overflow error occurred before sum across PEs.")

    call EFP_sum_across_PEs(asum(G%jsg:G%jeg), G%jeg-G%jsg+1)

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_i_mean overflow error occurred during sum across PEs.")

    do j=js,je
      i_mean(j) = EFP_to_real(asum(j+jdg_off)) / real(G%ieg-G%isg+1)
    enddo
  endif

  if (rescale /= 1.0) then ; do j=js,je ; i_mean(j) = rescale*i_mean(j) ; enddo ; endif

  deallocate(asum)

end subroutine global_i_mean

!> Determine the global mean of a field along rows of constant j, returning it
!! in a 1-d array using the local indexing. This uses reproducing sums.
subroutine global_j_mean(array, j_mean, G, mask, scale, tmp_scale, unscale)
  type(ocean_grid_type),            intent(inout) :: G     !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: array !< The variable to integrate in arbitrary units [a],
                                                           !! or arbitrary rescaled units [A ~> a] if unscale
                                                           !! is present
  real, dimension(SZI_(G)),         intent(out)   :: j_mean !<  Global mean of array along its j-axis [a] or [A ~> a]
  real, dimension(SZI_(G),SZJ_(G)), &
                          optional, intent(in)    :: mask  !< An array used for weighting the j-mean [nondim]
  real,                   optional, intent(in)    :: scale !< A rescaling factor for the output variable [a A-1 ~> 1]
                                                           !! that converts it back to unscaled (e.g., mks)
                                                           !! units to enable the use of the reproducing sums
  real,                   optional, intent(in)    :: tmp_scale !< A rescaling factor for the internal
                                                           !! calculations that is removed from the output [a A-1 ~> 1]
  real,                   optional, intent(in)    :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                           !! that converts it back to unscaled (e.g., mks)
                                                           !! units to enable the use of the reproducing sums.
                                                           !! Here scale and unscale are synonymous, but unscale
                                                           !! is preferred and takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  type(EFP_type), allocatable, dimension(:) :: asum      ! The masked sum of the variable in each row [a]
  type(EFP_type), allocatable, dimension(:) :: mask_sum  ! The sum of the mask values in each row [nondim]
  real :: mask_sum_r ! The sum of the mask values in a row [nondim]
  real :: scalefac   ! A scaling factor for the variable [a A-1 ~> 1]
  real :: rescale    ! A factor for redoing any internal rescaling before output [A a-1 ~> 1]
  integer :: is, ie, js, je, idg_off, jdg_off
  integer :: i, j

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  idg_off = G%idg_offset ; jdg_off = G%jdg_offset

  scalefac = 1.0
  if (present(unscale)) then ; scalefac = unscale
  elseif (present(scale)) then ; scalefac = scale ; endif

  rescale = 1.0
  if (present(tmp_scale)) then ; if (tmp_scale /= 0.0) then
    scalefac = scalefac * tmp_scale ; rescale = 1.0 / tmp_scale
  endif ; endif
  call reset_EFP_overflow_error()

  allocate(asum(G%isg:G%ieg))
  if (present(mask)) then
    allocate (mask_sum(G%isg:G%ieg))

    do i=G%isg,G%ieg
      asum(i) = real_to_EFP(0.0) ; mask_sum(i) = real_to_EFP(0.0)
    enddo

    do i=is,ie ; do j=js,je
      asum(i+idg_off) = asum(i+idg_off) + real_to_EFP(scalefac*array(i,j)*mask(i,j))
      mask_sum(i+idg_off) = mask_sum(i+idg_off) + real_to_EFP(mask(i,j))
    enddo ; enddo

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_j_mean overflow error occurred before sums across PEs.")

    call EFP_sum_across_PEs(asum(G%isg:G%ieg), G%ieg-G%isg+1)
    call EFP_sum_across_PEs(mask_sum(G%isg:G%ieg), G%ieg-G%isg+1)

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_j_mean overflow error occurred during sums across PEs.")

    do i=is,ie
      mask_sum_r = EFP_to_real(mask_sum(i+idg_off))
      if (mask_sum_r == 0.0 ) then ; j_mean(i) = 0.0 ; else
        j_mean(i) = EFP_to_real(asum(i+idg_off)) / mask_sum_r
      endif
    enddo

    deallocate(mask_sum)
  else
    do i=G%isg,G%ieg ; asum(i) = real_to_EFP(0.0) ; enddo

    do i=is,ie ; do j=js,je
      asum(i+idg_off) = asum(i+idg_off) + real_to_EFP(scalefac*array(i,j))
    enddo ; enddo

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_j_mean overflow error occurred before sum across PEs.")

    call EFP_sum_across_PEs(asum(G%isg:G%ieg), G%ieg-G%isg+1)

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_j_mean overflow error occurred during sum across PEs.")

    do i=is,ie
      j_mean(i) = EFP_to_real(asum(i+idg_off)) / real(G%jeg-G%jsg+1)
    enddo
  endif

  if (rescale /= 1.0) then ; do i=is,ie ; j_mean(i) = rescale*j_mean(i) ; enddo ; endif

  deallocate(asum)

end subroutine global_j_mean

!> Adjust 2d array such that area mean is zero without moving the zero contour
subroutine adjust_area_mean_to_zero(array, G, scaling, unit_scale, unscale)
  type(ocean_grid_type),            intent(in)    :: G       !< Grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: array   !< 2D array to be adjusted in  arbitrary units [a],
                                                             !! or arbitrary rescaled units [A ~> a] if unscale
                                                             !! is present
  real, optional,                   intent(out)   :: scaling !< The scaling factor used [nondim]
  real,                   optional, intent(in)    :: unit_scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                             !! that converts it back to unscaled (e.g., mks)
                                                             !! units to enable the use of the reproducing sums
  real,                   optional, intent(in)    :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                             !! that converts it back to unscaled (e.g., mks)
                                                             !! units to enable the use of the reproducing sums.
                                                             !! Here unit_scale and unscale are synonymous, but unscale
                                                             !! is preferred and takes precedence if both are present.
  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  real :: posVals(G%isc:G%iec, G%jsc:G%jec) ! The positive values in a cell or 0 [A ~> a]
  real :: negVals(G%isc:G%iec, G%jsc:G%jec) ! The negative values in a cell or 0 [A ~> a]
  real :: areaXposVals(G%isc:G%iec, G%jsc:G%jec) ! The cell area integral of the positive values [L2 A ~> m2 a]
  real :: areaXnegVals(G%isc:G%iec, G%jsc:G%jec) ! The cell area integral of the negative values [L2 A ~> m2 a]
  type(EFP_type), dimension(2) :: areaInt_EFP ! An EFP version integral of the values on the current PE [m2 a]
  real :: scalefac  ! A scaling factor for the variable [a A-1 ~> 1]
  real :: areaIntPosVals, areaIntNegVals ! The global area integral of the positive and negative values [m2 a]
  real :: posScale, negScale ! The scaling factor to apply to positive or negative values [nondim]
  integer :: i,j

  scalefac = 1.0
  if (present(unscale)) then ; scalefac = unscale
  elseif (present(unit_scale)) then ; scalefac = unit_scale ; endif

  ! areaXposVals(:,:) = 0.  ! This zeros out halo points.
  ! areaXnegVals(:,:) = 0.  ! This zeros out halo points.

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    posVals(i,j) = max(0., array(i,j))
    areaXposVals(i,j) = G%areaT(i,j) * posVals(i,j)
    negVals(i,j) = min(0., array(i,j))
    areaXnegVals(i,j) = G%areaT(i,j) * negVals(i,j)
  enddo ; enddo

  ! Combining the sums like this avoids separate blocking global sums.
  areaInt_EFP(1) = reproducing_sum_EFP( areaXposVals, only_on_PE=.true., unscale=scalefac*G%US%L_to_m**2 )
  areaInt_EFP(2) = reproducing_sum_EFP( areaXnegVals, only_on_PE=.true., unscale=scalefac*G%US%L_to_m**2 )
  call EFP_sum_across_PEs(areaInt_EFP, 2)
  areaIntPosVals = EFP_to_real( areaInt_EFP(1) )
  areaIntNegVals = EFP_to_real( areaInt_EFP(2) )

  posScale = 0.0 ; negScale = 0.0
  if ((areaIntPosVals>0.).and.(areaIntNegVals<0.)) then ! Only adjust if possible
    if (areaIntPosVals>-areaIntNegVals) then ! Scale down positive values
      posScale = - areaIntNegVals / areaIntPosVals
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        array(i,j) = (posScale * posVals(i,j)) + negVals(i,j)
      enddo ; enddo
    elseif (areaIntPosVals<-areaIntNegVals) then ! Scale down negative values
      negScale = - areaIntPosVals / areaIntNegVals
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        array(i,j) = posVals(i,j) + (negScale * negVals(i,j))
      enddo ; enddo
    endif
  endif
  if (present(scaling)) scaling = posScale - negScale

end subroutine adjust_area_mean_to_zero


!> Find the global maximum and minimum of a tracer array and return the locations of the extrema.
!! When there multiple cells with the same extreme values, the reported locations are from the
!! uppermost layer where they occur, and then from the logically northernmost and then eastermost
!! such location on the unrotated version of the grid within that layer.  Only ocean points (as
!! indicated by a positive value of G%mask2dT) are evaluated, and if there are no ocean points
!! anywhere in the domain, the reported extrema and their locations are all returned as 0.
subroutine array_global_min_max(tr_array, G, nk, g_min, g_max, &
                                xgmin, ygmin, zgmin, xgmax, ygmax, zgmax, unscale)
  integer,                      intent(in)  :: nk    !< The number of vertical levels
  type(ocean_grid_type),        intent(in)  :: G     !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G),nk), intent(in)  :: tr_array !< The tracer array to search for
                                                     !! extrema in arbitrary concentration units [CU ~> conc]
  real,                         intent(out) :: g_min !< The global minimum of tr_array, either in
                                                     !! the same units as tr_array [CU ~> conc] or in
                                                     !! unscaled units if unscale is present [conc]
  real,                         intent(out) :: g_max !< The global maximum of tr_array, either in
                                                     !! the same units as tr_array [CU ~> conc] or in
                                                     !! unscaled units if unscale is present [conc]
  real,               optional, intent(out) :: xgmin !< The x-position of the global minimum in the
                                                     !! units of G%geoLonT, often [degrees_E] or [km] or [m]
  real,               optional, intent(out) :: ygmin !< The y-position of the global minimum in the
                                                     !! units of G%geoLatT, often [degrees_N] or [km] or [m]
  real,               optional, intent(out) :: zgmin !< The z-position of the global minimum [layer]
  real,               optional, intent(out) :: xgmax !< The x-position of the global maximum in the
                                                     !! units of G%geoLonT, often [degrees_E] or [km] or [m]
  real,               optional, intent(out) :: ygmax !< The y-position of the global maximum in the
                                                     !! units of G%geoLatT, often [degrees_N] or [km] or [m]
  real,               optional, intent(out) :: zgmax !< The z-position of the global maximum [layer]
  real,               optional, intent(in)  :: unscale !< A factor to use to undo any scaling of
                                                     !! the input tracer array [conc CU-1 ~> 1]

  ! Local variables
  real    :: tmax, tmin      ! Maximum and minimum tracer values, in the same units as tr_array [CU ~> conc]
  integer :: ijk_min_max(2)  ! Integers encoding the global grid positions of the global minimum and maximum values
  real    :: xyz_min_max(6)  ! A single array with the x-, y- and z-positions of the minimum and
                             ! maximum values in units that vary between the array elements [various]
  logical :: valid_PE        ! True if there are any valid points on the local PE.
  logical :: find_location   ! If true, report the locations of the extrema
  integer :: ijk_loc_max     ! An integer encoding the global grid position of the maximum tracer value on this PE
  integer :: ijk_loc_min     ! An integer encoding the global grid position of the minimum tracer value on this PE
  integer :: ijk_loc_here    ! An integer encoding the global grid position of the current grid point
  integer :: itmax, jtmax, ktmax, itmin, jtmin, ktmin
  integer :: i, j, k, isc, iec, jsc, jec

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  find_location = (present(xgmin) .or. present(ygmin) .or. present(zgmin) .or. &
                   present(xgmax) .or. present(ygmax) .or. present(zgmax))

  ! The initial values set here are never used if there are any valid points.
  tmax = -huge(tmax) ; tmin = huge(tmin)

  if (find_location) then
    ! Find the maximum and minimum tracer values on this PE and their locations.
    valid_PE = .false.
    itmax = 0 ; jtmax = 0 ; ktmax = 0 ; ijk_loc_max = 0
    itmin = 0 ; jtmin = 0 ; ktmin = 0 ; ijk_loc_min = 0
    do k=1,nk ; do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j) > 0.0) then
      valid_PE = .true.
      if (tr_array(i,j,k) > tmax) then
        tmax = tr_array(i,j,k)
        itmax = i ; jtmax = j ; ktmax = k
        ijk_loc_max = ijk_loc(i, j, k, nk, G%HI)
      elseif ((tr_array(i,j,k) == tmax) .and. (k <= ktmax)) then
        ijk_loc_here = ijk_loc(i, j, k, nk, G%HI)
        if (ijk_loc_here > ijk_loc_max) then
          itmax = i ; jtmax = j ; ktmax = k
          ijk_loc_max = ijk_loc_here
        endif
      endif
      if (tr_array(i,j,k) < tmin) then
        tmin = tr_array(i,j,k)
        itmin = i ; jtmin = j ; ktmin = k
        ijk_loc_min = ijk_loc(i, j, k, nk, G%HI)
      elseif ((tr_array(i,j,k) == tmin) .and. (k <= ktmin)) then
        ijk_loc_here = ijk_loc(i, j, k, nk, G%HI)
        if (ijk_loc_here > ijk_loc_min) then
          itmin = i ; jtmin = j ; ktmin = k
          ijk_loc_min = ijk_loc_here
        endif
      endif
    endif ; enddo ; enddo ; enddo
  else
    ! Only the maximum and minimum values are needed, and not their positions.
    do k=1,nk ; do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j) > 0.0) then
      if (tr_array(i,j,k) > tmax) tmax = tr_array(i,j,k)
      if (tr_array(i,j,k) < tmin) tmin = tr_array(i,j,k)
    endif ; enddo ; enddo ; enddo
  endif

  ! Find the global maximum and minimum tracer values.
  g_max = tmax ; g_min = tmin
  call max_across_PEs(g_max)
  call min_across_PEs(g_min)

  if (find_location) then
    if (g_max < g_min) then
      ! This only occurs if there are no unmasked points anywhere in the domain.
      xyz_min_max(:) = 0.0
    else
      ! Find the global indices of the maximum and minimum locations.  This can
      ! occur on multiple PEs.
      ijk_min_max(1:2) = 0
      if (valid_PE) then
        if (g_min == tmin) ijk_min_max(1) = ijk_loc_min
        if (g_max == tmax) ijk_min_max(2) = ijk_loc_max
      endif
      ! If MOM6 supported taking maxima on arrays of integers, these could be combined as:
      ! call max_across_PEs(ijk_min_max, 2)
      call max_across_PEs(ijk_min_max(1))
      call max_across_PEs(ijk_min_max(2))

      ! Set the positions of the extrema if they occur on this PE.  This will only
      ! occur on a single PE.
      xyz_min_max(1:6) = -huge(xyz_min_max)  ! These huge negative values are never selected by max_across_PEs.
      if (valid_PE) then
        if (ijk_min_max(1) == ijk_loc_min) then
          xyz_min_max(1) = G%geoLonT(itmin,jtmin)
          xyz_min_max(2) = G%geoLatT(itmin,jtmin)
          xyz_min_max(3) = real(ktmin)
        endif
        if (ijk_min_max(2) == ijk_loc_max) then
          xyz_min_max(4) = G%geoLonT(itmax,jtmax)
          xyz_min_max(5) = G%geoLatT(itmax,jtmax)
          xyz_min_max(6) = real(ktmax)
        endif
      endif

      call max_across_PEs(xyz_min_max, 6)
    endif

    if (present(xgmin)) xgmin = xyz_min_max(1)
    if (present(ygmin)) ygmin = xyz_min_max(2)
    if (present(zgmin)) zgmin = xyz_min_max(3)
    if (present(xgmax)) xgmax = xyz_min_max(4)
    if (present(ygmax)) ygmax = xyz_min_max(5)
    if (present(zgmax)) zgmax = xyz_min_max(6)
  endif

  if (g_max < g_min) then
    ! There are no unmasked points anywhere in the domain.
    g_max = 0.0 ; g_min = 0.0
  endif

  if (present(unscale)) then
    ! Rescale g_min and g_max, perhaps changing their units from [CU ~> conc] to [conc]
    g_max = unscale * g_max
    g_min = unscale * g_min
  endif

end subroutine array_global_min_max

! Return a positive integer encoding the rotationally invariant global position of a tracer cell
function ijk_loc(i, j, k, nk, HI)
  integer,              intent(in) :: i   !< Local i-index
  integer,              intent(in) :: j   !< Local j-index
  integer,              intent(in) :: k   !< Local k-index
  integer,              intent(in) :: nk  !< Range of k-index, used to pick out a low-k position.
  type(hor_index_type), intent(in) :: HI  !< Horizontal index ranges
  integer :: ijk_loc  ! An integer encoding the cell position in the global grid.

  ! Local variables
  integer :: ig, jg  ! Global index values with a global computational domain start value of 1.
  integer :: ij_loc  ! The encoding of the horizontal position
  integer :: qturns  ! The number of counter-clockwise quarter turns of the grid that have to be undone

  ! These global i-grid positions run from 1 to HI%niglobal, and analogously for jg.
  ig = i + HI%idg_offset + (1 - HI%isg)
  jg = j + HI%jdg_offset + (1 - HI%jsg)

  ! Compensate for the rotation of the model grid to give a rotationally invariant encoding.
  qturns = modulo(HI%turns, 4)
  if (qturns == 0) then
    ij_loc = ig + HI%niglobal * jg
  elseif (qturns == 1) then
    ij_loc = jg + HI%njglobal * ((HI%niglobal+1)-ig)
  elseif (qturns == 2) then
    ij_loc = ((HI%niglobal+1)-ig) + HI%niglobal * ((HI%njglobal+1)-jg)
  elseif (qturns == 3) then
    ij_loc = ((HI%njglobal+1)-jg) + HI%njglobal * ig
  endif

  ijk_loc = ij_loc + (HI%niglobal*HI%njglobal) * (nk-k)

end function ijk_loc


end module MOM_spatial_means
